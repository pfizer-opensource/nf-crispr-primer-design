"""Additional output to support NGS."""
import logging
from pathlib import Path

import yaml

from .exceptions import PrimerDesignNGSPrimerError
from .subcommands import bed_intersect_cmd
from .primers import get_amplicon_name, make_primers_bed_table


_logger = logging.getLogger(__name__)

def get_amplicon_coding_coordinates(primers, ref_files, config):
    """Get coordinates for overlap between amplicons and coding regions."""
    _logger.info('Finding overlap between amplicons and coding exons...')
    coding_bed = ref_files['coding_bed'][config['coding']]
    amplicon_bed = '\n'.join(
        '\t'.join(row) for row in make_primers_bed_table(primers, header=False)
    )
    if not amplicon_bed:
        return
    amplicon_exons = bed_intersect_cmd(coding_bed, config)
    for region in amplicon_exons.communicate(f'{amplicon_bed}\n'):
        fields = region.split('\t')
        name_parts = fields[3].split('-pair')
        id = name_parts[0]
        pair = 0 if len(name_parts) == 1 else int(name_parts[1]) - 1
        start = int(fields[1])
        end = int(fields[2])
        # Store start/end coordinates in a set because some genes have multiple
        # overlapping transcripts with the same exons and and we only want to
        # keep one copy of the same coding region.
        try:
            primers[id]['pairs'][pair]['coding'].add((start, end))
        except KeyError:
            primers[id]['pairs'][pair]['coding'] = {(start, end)}

def get_amplicons_path(amplicons_file, output_dir):
    amplicons_dir = Path(amplicons_file).parent
    if amplicons_dir == Path():
        amplicons_dir = output_dir
    return amplicons_dir / Path(amplicons_file).name

def make_amplicons_yaml(primers, ref_files, config):
    """Return an amplicons YAML document for the NGS CRISPR Pipeline."""
    if config['coding'] is not None:
        get_amplicon_coding_coordinates(primers, ref_files, config)
    _logger.info('Generating amplicons YAML...')
    amplicons = {}
    multiple_guide_pairs = False
    # Build YAML dict of all amplicons
    for id, guide_primers in primers.items():
        for pair in guide_primers['pairs']:
            # Genomic sequence of amplicon
            seq = ''.join((guide_primers['guide']['L_flank_seq'],
                           guide_primers['guide']['genomic_seq'],
                           guide_primers['guide']['R_flank_seq']))
            amplicon_seq = seq[pair['F_primer_offset']:pair['R_primer_offset'] + 1]
            amplicon_name = get_amplicon_name(id, pair)
            if pair['total_pairs'] > 1:
                multiple_guide_pairs = True
            # Genomic sequence of coding regions
            coding_seqs = []
            try:
                for start, end in pair['coding']:
                    start_offset = start - pair['F_primer_start']
                    end_offset = end - pair['F_primer_start']
                    coding_seqs.append(amplicon_seq[start_offset:end_offset])
            except KeyError:
                pass
            amplicons[amplicon_name] = {
                'seq': amplicon_seq,
                'guide': guide_primers['guide']['seq'],
            }
            if coding_seqs:
                amplicons[amplicon_name]['coding'] = ','.join(coding_seqs)
    if multiple_guide_pairs:
        _logger.warning('Multiple amplicons designed per target site. '
                        'Amplicons YAML file may need to be edited to only '
                        'include the desired primer pairs.')
    return yaml.dump(amplicons, Dumper=yaml.SafeDumper, sort_keys=False)

def get_ngs_adapters(config):
    """Return tuple of adapter name and NGS adapters dict from config."""
    adapter_name = config['ngs_adapters']
    try:
        return (adapter_name, config['ngs']['adapters'][adapter_name])
    except KeyError:
        raise PrimerDesignNGSPrimerError(
            f'Invalid name of NGS adapters: {adapter_name}'
        ) from None

def make_ngs_primers_table(primers, config, input=None):
    """Return a generator that provides each row of the NGS primer table as a
    list of fields."""
    adapter_name, adapters = get_ngs_adapters(config)
    header_names = ['Amplicon name']
    if input is None:
        input_rows = None
        header_names.append('Guide sequence')
    else:
        input_rows = iter(input)
        header_names.extend(next(input_rows).rstrip('\r\n').split('\t'))
    header_names.extend((
        'Target-specific F primer',
        'Target-specific R primer',
        f'{adapter_name} pair name',
        f'{adapter_name} F primer sequence',
        f'{adapter_name} R primer sequence',
        'F primer length',
        'R primer length',
        'Amplicon length',
        'Flag duplicates',
    ))
    yield header_names
    for id, guide_primers in primers.items():
        if input_rows is None:
            guide_info = (guide_primers['guide']['seq'], )
        else:
            guide_info = next(input_rows).rstrip('\r\n').split('\t')
        if guide_primers['pairs']:
            for pair in guide_primers['pairs']:
                amplicon_name = get_amplicon_name(id, pair)
                # Primer pair name is the amplicon name with the optional adapter suffix
                if adapters['suffix']:
                    pair_name = f'{amplicon_name} {adapters["suffix"]}'
                else:
                    pair_name = amplicon_name
                # Calculate number of N bases to insert between the adapter and the
                # target-specific sequence for forward and reverse primers
                forward_num_n = max(config['ngs']['min_n'],
                                    min(config['ngs']['max_n'],
                                        config['ngs']['max_length']
                                        - len(adapters['forward'])
                                        - len(pair['F_primer_seq'])))
                reverse_num_n = max(config['ngs']['min_n'],
                                    min(config['ngs']['max_n'],
                                        config['ngs']['max_length']
                                        - len(adapters['reverse'])
                                        - len(pair['R_primer_seq'])))
                # NGS primers are built from the adapter, an optional run of Ns,
                # and the target-specific sequence for each direction
                forward_seq = ''.join((adapters['forward'],
                                       'N' * forward_num_n,
                                       pair['F_primer_seq']))
                reverse_seq = ''.join((adapters['reverse'],
                                       'N' * reverse_num_n,
                                       pair['R_primer_seq']))
                amplicon_length = (pair['amplicon_size']
                                   - len(pair['F_primer_seq'])
                                   - len(pair['R_primer_seq'])
                                   + len(forward_seq)
                                   + len(reverse_seq))
                # Yield row of table for this amplicon
                row = [amplicon_name]
                row.extend(guide_info)
                row.extend((
                    pair['F_primer_seq'],
                    pair['R_primer_seq'],
                    pair_name,
                    forward_seq.upper(),
                    reverse_seq.upper(),
                    len(forward_seq),
                    len(reverse_seq),
                    amplicon_length,
                    pair['duplicates'],
                ))
                yield row
        else:
            # No primer pairs for this row, so output guide ID (as Amplicon name) and
            # guide info with empty fields for remaining columns
            row = [id]
            row.extend(guide_info)
            row.extend([''] * (len(header_names) - len(guide_info) - 1))
            yield row
