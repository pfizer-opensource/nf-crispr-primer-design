from pathlib import Path
from pprint import pformat
import os
import sys

from confuse import ConfigError, ConfigSource, RootView
from confuse.templates import MappingValues, Optional, StrSeq
import yaml


# Subcommand configuration template
subcommand_template = {
    'cmd': StrSeq(), # List of strings, required
    'lmod': Optional(StrSeq()), # List of strings, optional, default to None
}

# Primer design configuration template
pd_config_template = {
    'amplicons': [bool, str, Path], # Boolean, string or path, required
    'coding': Optional(str), # String, optional, default to None
    'ngs_adapters': Optional(str), # String, optional, default to None
    'include_input': bool, # Boolean, required
    'input_format': {
        'id_col': int, # Integer, required
        'seq_col': int, # Integer, required
        'chr_col': Optional(int), # Integer, optional, default to None
        'pos_col': Optional(int), # Integer, optional, default to None
        'strand_col': Optional(int) # Integer, optional, default to None
    },
    'parameters': {
        'size_range': StrSeq(split=True), # List of strings or whitespace separate string, required
        'size_opt': int, # Integer, required
        'target': int, # Integer, required
        'flanking': int, # Integer, required
        'num_pairs': int, # Integer, required
        'primer_opt': int, # Integer, required
        'primer_min': int, # Integer, required
        'primer_max': int, # Integer, required
        'max_poly_x': int # Integer, required
    },
    'ngs': {
        'adapters': MappingValues({
            'forward': str, # String, required
            'reverse': str, # String, required
            'suffix': '', # String, optional, default to empty string
        }),
        'min_n': 0, # Integer, optional, default to 0
        'max_n': 0, # Integer, optional, default to 0
        'max_length': Optional(1000), # Integer, optional, default to 1000 (ie, any length)
    },
    'references': MappingValues({
        'fasta': [str, Path], # String or path, required
        'coding_bed': MappingValues(Optional(str)), # Mapping sequence of strings
        'genome_base': Optional(str), # String, optional, default to None
        'chrom_sizes': Optional(str), # String, optional, default to None
        'priority': 0 # Integer, optional, default to 0
    }),
    'subcommands': {
        'bowtie2': subcommand_template,
        'samtools': subcommand_template,
        'bedtools': subcommand_template,
    },
}


# Built-in default config
builtin_config = {
    'amplicons': False,
    'include_input': False,
    'input_format': {
        'id_col': 1,
        'seq_col': 2
    },
    'parameters': {
        'size_range': ['125-280'],
        'size_opt': 250,
        'target': 100,
        'flanking': 500,
        'num_pairs': 1,
        'primer_opt': 20,
        'primer_min': 18,
        'primer_max': 24,
        'max_poly_x': 4
    },
    'ngs': {},
    'references': {},
    'subcommands': {
        'bowtie2': {'cmd': 'bowtie2'},
        'samtools': {'cmd': 'samtools'},
        'bedtools': {'cmd': 'bedtools'},
    },
}


class ConfigSourceNamed(ConfigSource):
    """Extend ConfigSource with source name and error message."""
    def __init__(self, value, filename=None, default=False,
                 name=None, error=None):
        super().__init__(value, filename=filename, default=default)
        self.name = name
        self.error = error


def pd_yaml_config(path, name=None):
    """Produce a configuration source from a primer design yaml file."""
    if path is None:
        return ConfigSourceNamed({}, name=name, error='Not provided')
    
    try:
        return ConfigSourceNamed(load_pd_yaml(path), str(path), name=name)
    except OSError:
        error = 'Not present'
    except (yaml.YAMLError, ConfigError, ValueError):
        error = 'Error parsing configuration file'
        print('Warning: {0}: {1}'.format(error, path), file=sys.stderr)

    return ConfigSourceNamed({}, str(path), name=name, error=error)


def load_pd_yaml(path):
    """Load a primer design yaml file.
    If yaml contains a "primer_design" section, use that as the configuration.
    Otherwise use the entire yaml document as the configuration.
    """
    with open(path, 'r') as yaml_file:
        yaml_cfg = yaml.safe_load(yaml_file)
        if 'primer_design' in yaml_cfg:
            yaml_cfg = yaml_cfg['primer_design']
        return yaml_cfg


def load_default_config():
    """Load default configuration from multiple cascading sources.
    Default configuration will be loaded from the following sources in
    descending order of precedence:
      - A file specified in `$PRIMER_DESIGN_CONFIG`.
      - A file `primer_design.yaml` in the current working directory.
      - A file `$XDG_CONFIG_HOME/CRISPR_primer_design/config.yaml`,
        which defaults to `~/.config/CRISPR_primer_design/config.yaml`.
      - A file `config.yaml` in the installation directory (next to the script).
      - Hard-coded built-in default configuration.
    """
    xdg_config_home = Path(os.environ.get('XDG_CONFIG_HOME', '~/.config'))
     
    yaml_sources = [
        ('$PRIMER_DESIGN_CONFIG', os.environ.get('PRIMER_DESIGN_CONFIG')),
        ('working directory', Path.cwd() / 'primer_design.yaml'),
        ('$XDG_CONFIG_HOME', xdg_config_home.expanduser() / 'CRISPR_primer_design/config.yaml'),
        ('installation directory', Path(sys.argv[0]).resolve().parent / 'config.yaml'),
    ]
    
    config = RootView([])
    for name, yaml_path in yaml_sources:
        config.add(pd_yaml_config(yaml_path, name))
    
    config.add(ConfigSourceNamed(builtin_config, name='built-in configuration'))
    
    return config


def print_config(config, template=None, sources=True, final=True, file=sys.stdout):
    """Print configuration sources and/or final config.
    Pretty print a set of cascading configuration sources and the final
    configuration as validated by an optional template.
    """
    if sources:
        print('# Configuration sources (descending precedence):', file=file)
        for source in config.sources:
            name = getattr(source, 'name', '')
            error = getattr(source, 'error', None)
            print('## {0}{1}\n{2}\n'.format(
                name+' ' if name else '',
                '('+source.filename+')' if source.filename else '',
                error or pformat(source)
            ), file=file)
    
    if sources and final:
        print('', file=file)
    
    if final:
        print('# Final configuration:\n{0}'.format(
            pformat(config.get(template))
        ), file=file)
