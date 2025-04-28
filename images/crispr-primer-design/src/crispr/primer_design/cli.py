import argparse
import logging
from pathlib import Path
import re
import sys

from yaml import YAMLError

from .config import ConfigError, ConfigSourceNamed, load_default_config, load_pd_yaml, pd_config_template, print_config
from .exceptions import PrimerDesignError, PrimerDesignIOError
from .guides import prepare_guides
from .io import InputReader, prep_output_dir, write_tsv_table
from .ngs import get_amplicons_path, get_ngs_adapters, make_amplicons_yaml, make_ngs_primers_table
from .primers import design_primers, designer, make_primers_table, make_primers_bed_table
from .ref_genomes import get_reference_files, get_ref_names


__version__ = '0.3.11'
program_name = 'CRISPR primer design'
program_info = f'{program_name} v{__version__}, using {designer}'
_logger = logging.getLogger(__name__)

def handle_args(defaults):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=program_info,
                                     allow_abbrev=False, add_help=False)
    
    # Required parameters
    req_args = parser.add_argument_group('Required parameters')
    req_args.add_argument('-i', '--input', type=str,
                          metavar='FILE', required=True,
                          help='''Input file (tab-delimited text) with guide ID
                               and sequence columns. Use "-" to read input from
                               STDIN.''')
    req_args.add_argument('-o', '--output', type=Path,
                          metavar='FILE', required=True,
                          help='Output file for primer design.')
    req_args.add_argument('-r', '--reference', type=str,
                          metavar='NAME', required=True,
                          help='''Name of reference genome. Use --list-refs
                               to see options.''')
    
    # Optional parameters
    opt_args = parser.add_argument_group('Optional parameters')
    opt_args.add_argument('-a', '--amplicons', type=Path, metavar='FILE',
                          nargs='?', const=True, default=None,
                          help='Name of optional amplicons YAML file to produce.')
    opt_args.add_argument('--no-amplicons', dest='amplicons', action='store_false',
                          help=f'''Do not produce an amplicons YAML file.
                               {"(default)" if not defaults.amplicons else ""}''')
    opt_args.add_argument('-c', '--coding', type=str, metavar='NAME',
                          help='''Name of coding exon definitions to use when
                               making amplicons file. Use --list-refs to see
                               options.''')
    opt_args.add_argument('--ngs-adapters', type=str, metavar='NAME',
                          help='''Name of NGS adapters to append to primers.
                               Use --list-adapters to see options.''')
    opt_args.add_argument('--include-input', action='store_const', const=True,
                          help=f'''Include all columns from the input file in
                               the output file.
                               {"(default)" if defaults.include_input else ""}''')
    opt_args.add_argument('--no-include-input', action='store_false',
                          dest='include_input',
                          help=f'''Do not include the additional input columns
                               in the output.
                               {"" if defaults.include_input else "(default)"}''')
    opt_args.add_argument('--id-col', type=int, dest='input_format.id_col',
                          metavar='NUM',
                          help=f'''Column of input containing the guide ID.
                               (default: {defaults.input_format.id_col})''')
    opt_args.add_argument('--seq-col', type=int, dest='input_format.seq_col',
                          metavar='NUM',
                          help=f'''Column of input containing the guide sequence.
                               (default: {defaults.input_format.seq_col})''')
    opt_args.add_argument('--chr-col', type=int, dest='input_format.chr_col',
                          metavar='NUM',
                          help=f'''Column of input containing the chromosome of
                               the guide if already determined.
                               (default: {defaults.input_format.chr_col})''')
    opt_args.add_argument('--pos-col', type=int, dest='input_format.pos_col',
                          metavar='NUM',
                          help=f'''Column of input containing the genomic
                               position of the guide if already determined.
                               (default: {defaults.input_format.pos_col})''')
    opt_args.add_argument('--strand-col', type=int, dest='input_format.strand_col',
                          metavar='NUM',
                          help=f'''Column of input containing the genomic
                               strand of the guide if already determined.
                               (default: {defaults.input_format.strand_col})''')
    opt_args.add_argument('--size-range', type=size_range, dest='parameters.size_range',
                          metavar='MIN-MAX', action='append',
                          help=f'''Range of amplicon sizes ('min-max') to use
                               for primer design. This option can be used
                               multiple times, with ranges listed earlier
                               preferred over later ranges.
                               (default: {', '.join(defaults.parameters.size_range)})''')
    opt_args.add_argument('--size-opt', type=int, dest='parameters.size_opt',
                          metavar='NUM',
                          help=f'''Optimal size for the amplicon.
                               (default: {defaults.parameters.size_opt})''')
    opt_args.add_argument('--target', type=int, dest='parameters.target',
                          metavar='NUM',
                          help=f'''Number of base pairs centered on the guide
                               that must be included within the amplicon.
                               (default: {defaults.parameters.target})''')
    opt_args.add_argument('--flanking', type=int, dest='parameters.flanking',
                          metavar='NUM',
                          help=f'''Number of base pairs on either side of guide
                               to search for primer designs.
                               (default: {defaults.parameters.flanking})''')
    opt_args.add_argument('--num-pairs', type=int, dest='parameters.num_pairs',
                          metavar='NUM',
                          help=f'''Number of primer pairs to design per guide.
                               (default: {defaults.parameters.num_pairs})''')
    opt_args.add_argument('--primer-opt', type=int, dest='parameters.primer_opt',
                          metavar='NUM',
                          help=f'''Optimum primer length.
                               (default: {defaults.parameters.primer_opt})''')
    opt_args.add_argument('--primer-min', type=int, dest='parameters.primer_min',
                          metavar='NUM',
                          help=f'''Minimum primer length.
                               (default: {defaults.parameters.primer_min})''')
    opt_args.add_argument('--primer-max', type=int, dest='parameters.primer_max',
                          metavar='NUM',
                          help=f'''Maximum primer length.
                               (default: {defaults.parameters.primer_max})''')
    opt_args.add_argument('--max-poly-x', type=int, dest='parameters.max_poly_x',
                          metavar='NUM',
                          help=f'''Maximum allowed length of a mononucleotide
                               repeat in a primer.
                               (default: {defaults.parameters.max_poly_x})''')
    opt_args.add_argument('--ngs-min-n', type=int, dest='ngs.min_n',
                          metavar='NUM',
                          help=f'''Minimum number of N bases to include in NGS
                               primer design to add early read complexity.
                               (default: {defaults.ngs.min_n})''')
    opt_args.add_argument('--ngs-max-n', type=int, dest='ngs.max_n',
                          metavar='NUM',
                          help=f'''Maximum number of N bases to include in NGS
                               primer design to add early read complexity.
                               (default: {defaults.ngs.min_n})''')
    opt_args.add_argument('--ngs-max-length', type=int, dest='ngs.max_length',
                          metavar='NUM',
                          help=f'''Maximum NGS primer length with adapter and
                               N's for early read complexity. --ngs-min-n
                               option takes precedence.
                               (default: {defaults.ngs.max_length})''')
    
    # Other options
    other = parser.add_argument_group('Other options')
    other.add_argument('--config', type=pd_yaml_config_arg, metavar='FILE',
                       help='Use configuration file in YAML format.')
    other.add_argument('--list-refs', action=ListRefs,
                       help='List available reference and coding exon names.')
    other.add_argument('--list-adapters', action=ListAdapters,
                       help='List available NGS adapters.')
    other.add_argument('-p', '--print-config', action=PrintConfig,
                       help='Print configuration.')
    other.add_argument('-d', '--debug', action='store_true', dest='debug',
                       help='Enable debug mode.')
    other.add_argument('-V', '--version', action='version', version=__version__,
                       help='Show program version number and exit.')
    other.add_argument('-h', '--help', action='help',
                       help='Show this help message and exit.')
    
    return parser.parse_args()


def size_range(size_range):
    """Validate that the size-range option matches the expected min-max pattern."""
    m = re.match(r'\s*((\d+)-(\d+))\s*$', size_range)
    if m and int(m.group(3)) > int(m.group(2)):
        return m.group(1)
    else:
        raise(argparse.ArgumentTypeError('Invalid value does not match format "min-max": {0}'.format(size_range)))


def pd_yaml_config_arg(path):
    """Validate the primer design yaml file given as --config argument.
    Returns a configuration source or raises an argparse error.
    """
    try:
        return ConfigSourceNamed(load_pd_yaml(path), str(path),
                                 name='--config argument file')
    except OSError:
        raise argparse.ArgumentTypeError('Could not read file {0}'.format(path))
    except (YAMLError, ConfigError, ValueError):
        raise argparse.ArgumentTypeError('Error parsing configuration file {0}'.format(path))


class ListRefs(argparse.Action):
    """Action to list references in configuration."""
    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest,
                         default=default, nargs=0, help=help)
    
    def __call__(self, parser, namespace, values, option_string=None):
        config = load_default_config()
        if namespace.config:
            config.set(namespace.config)
        pd_config = validate_pd_config(config)
        ref_names = get_ref_names(pd_config.references)
        print('Available reference genomes and coding exon definitions')
        if ref_names:
            for name, coding in ref_names.items():
                coding_aliases = ', '.join([x for x in coding if x != 'none'])
                print('{0}: {1}'.format(name, coding_aliases or '(no coding exons available)'))
        else:
            print('No references defined in configuration')
        parser.exit()


class ListAdapters(argparse.Action):
    """Action to list NGS adapters in configuration."""
    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest,
                         default=default, nargs=0, help=help)
    
    def __call__(self, parser, namespace, values, option_string=None):
        config = load_default_config()
        if namespace.config:
            config.set(namespace.config)
        pd_config = validate_pd_config(config)
        adapter_names = pd_config.ngs.adapters.keys()
        if adapter_names:
            print('Available NGS adapter definitions:')
            for name in adapter_names:
                print(f'  - {name}')
        else:
            print('No NGS adapters defined in configuration.')
        parser.exit()


class PrintConfig(argparse.Action):
    """Action to print current configuration."""
    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest,
                         default=default, nargs=0, help=help)
    
    def __call__(self, parser, namespace, values, option_string=None):
        config = load_default_config()
        if namespace.config:
            config.set(namespace.config)
            namespace.config = None
        config.set_args(namespace, dots=True)
        config.sources[0].name = 'command-line arguments'
        try:
            print_config(config, pd_config_template)
        except ConfigError as e:
            print('ERROR: Could not print final configuration due to validation error: {0}'.format(e),
                  file=sys.stderr)
        parser.exit()


def validate_pd_config(config):
    """Validate primer design configuration against template.
    If the configuration does not validate, exit with an error, but if the
    '-p' or '--print-config' options were passed, then print the sources.
    """
    try:
        return config.get(pd_config_template)
    except ConfigError as e:
        if '-p' in sys.argv or '--print-config' in sys.argv:
            print_config(config, pd_config_template, final=False)
            msg = 'Could not print final configuration due to validation error'
        else:
            msg = 'invalid configuration'
        _logger.error('{0}: {1}'.format(msg, e))
        sys.exit(1)


def main():
    """Primer design CLI entry point."""
    root_logger = logging.getLogger()
    log_handler = logging.StreamHandler(sys.stdout)
    log_handler.setFormatter(
        logging.Formatter('%(levelname)s: %(message)s')
    )
    root_logger.addHandler(log_handler)
    root_logger.setLevel(logging.INFO)

    config = load_default_config()
    defaults = validate_pd_config(config)
    
    args = handle_args(defaults)
    print(program_info, file=sys.stdout)
    if args.debug:
        log_handler.setFormatter(
            logging.Formatter('%(levelname)s:%(name)s: %(message)s')
        )
        root_logger.setLevel(logging.DEBUG)
        _logger.debug(f'config defaults:\n{defaults}')
    _logger.debug(f'arguments:\n{vars(args)}')
    if args.config:
        config.set(args.config)
        args.config = None
    config.set_args(args, dots=True)
    config.sources[0].name = 'command-line arguments'
    
    pd_config = validate_pd_config(config)
    _logger.debug(f'validated config:\n{pd_config}')
    
    try:
        input = InputReader(args.input, 'Enter guide sequences below:')
        output_dir = prep_output_dir(args.output)
        reference_files = get_reference_files(pd_config.references,
                                              args.reference,
                                              pd_config.coding)
        if pd_config.ngs_adapters:
            # Will raise an error early if ngs adapter name is not valid
            get_ngs_adapters(pd_config)
        guides = prepare_guides(input, reference_files, pd_config)
        primers = design_primers(guides, pd_config.parameters)
        _logger.info('Generating primer output...')
        write_tsv_table(args.output,
                        make_primers_table(primers,
                                           pd_config.parameters['num_pairs'],
                                           input if pd_config.include_input else None))
        _logger.info('Generating BED output...')
        bed_output = args.output.with_suffix('.bed')
        write_tsv_table(bed_output,
                           make_primers_bed_table(primers))
        if pd_config.ngs_adapters:
            ngs_filename = args.output.stem + '_NGS.txt'
            ngs_output = args.output.with_name(ngs_filename)
            _logger.info('Generating NGS primer output...')
            write_tsv_table(ngs_output,
                            make_ngs_primers_table(primers, pd_config,
                                                   input if pd_config.include_input else None))
        if pd_config.amplicons:
            amplicons_file = pd_config.amplicons if pd_config.amplicons is not True else 'amplicons.yaml'
            try:
                with open(get_amplicons_path(amplicons_file, output_dir), 'w') as yaml_out:
                    print(make_amplicons_yaml(primers, reference_files, pd_config),
                          end='', file=yaml_out)
            except OSError as e:
                raise PrimerDesignIOError(
                    f'Could not open amplicons output file "{amplicons_file}": {e}'
                ) from None
            
    except PrimerDesignError as e:
        _logger.error(str(e))
        sys.exit(1)
    else:
        _logger.info('Done!')
        sys.exit(0)
