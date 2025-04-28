"""Subcommand handling."""
import logging
import os
import subprocess
import sys


_logger = logging.getLogger(__name__)


class Subcommand:
    """A subcommand that processs data and communicates through stdin/stdout.
    """
    
    def __init__(self, command, *base_args, shell=False, text=True):
        self.command = command
        self.base_args = base_args
        self.shell = shell
        self.text = text
        self.returncode = None
        self._process = None
    
    def start(self, *args):
        """Start the subcommand."""
        if self._process is not None:
            _logger.warning('Subcommand already running: %s', self.command)
            return
        _args = [self.command]
        _args.extend(self.base_args)
        if args:
            _args.extend(args)
        self._process = subprocess.Popen(_args, stdin=subprocess.PIPE,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         shell=self.shell, text=self.text)
        self.returncode = None
    
    def stop(self):
        """Stop the subcommand, returning data from stdout and logging stderr."""
        if self._process is not None:
            stdout, stderr = self._process.communicate()
            self.returncode = self._process.returncode
            self._process = None
            if stderr:
                _logger.log(logging.WARNING if self.returncode else logging.INFO,
                            stderr.rstrip('\r\n'))
            return stdout
    
    def run(self, *args):
        """Run the subcommand and return any data."""
        self.start(*args)
        return self.stop()
    
    def communicate(self, input, args=(), end_of_record=None):
        """Send data to the process and return an output generator.
        The generator will yield each line of the process's output.
        The subcommand will be started with the provided `args` if it is not
        already running.
        If the `end_of_record` parameter is set, then the output will stop when
        this value is returned by the process, but the process will be kept
        running to accept more input. Otherwise, all of the available output
        will be returned and the process will be stopped.
        """
        if self._process is None:
            self.start(*args)
        _logger.debug(f'Input to subcommand:\n{input}')
        self._process.stdin.write(input)
        self._process.stdin.flush()
        if end_of_record is not None:
            output = self._process.stdout
        else:
            output = self.stop().splitlines(True)
        for line in output:
            if line == end_of_record:
                _logger.debug(f'Output from subcommand reached end of record')
                return
            _logger.debug(f'Output line from subcommand: {line}')
            yield line
        _logger.debug(f'Output from subcommand complete')


def load_modules(lmods):
    try:
        lmod_cmd = os.environ['LMOD_CMD']
    except KeyError:
        _logger.error('lmod command not found')
    lmod = Subcommand(lmod_cmd, 'python')
    exec(lmod.run('purge'))
    try:
        exec(lmod.run('load', *lmods))
    except NameError:
        pass

def guide_bed_cmd(reference_files, config):
    """Return a subcommand pipeline that takes in guide sequences and returns
    bed-formatted genomic coordinates.
    """
    commands = config['subcommands']
    modules = []
    for cmd in ('bowtie2', 'samtools', 'bedtools'):
        if commands[cmd]['lmod']:
            modules.extend(commands[cmd]['lmod'])
    if modules:
        _logger.debug(f'LMOD modules for guide to bed subcommand: {modules}')
        load_modules(modules)
    guide_align_cmd = f'{" ".join(commands["bowtie2"]["cmd"])} -x "{reference_files["genome_base"]}" -f -U -'
    sam2bam_cmd = f'{" ".join(commands["samtools"]["cmd"])} view -b -'
    bam2bed_cmd = f'{" ".join(commands["bedtools"]["cmd"])} bamtobed -i stdin'
    guide_bed_cmd = f'set -o pipefail; {guide_align_cmd} | {sam2bam_cmd} | {bam2bed_cmd}'
    _logger.debug(f'Guide to bed subcommand: {guide_bed_cmd}')
    return Subcommand(guide_bed_cmd, shell=True)

def guide_flank_cmd(reference_files, config):
    """Return a subcommand pipeline that takes in bed-formatted genomic
    coordinates and returns the flanking genomic DNA sequence.
    """
    if config['subcommands']['bedtools']['lmod']:
        _logger.debug('LMOD modules for guide flank subcommand: {}'.format(
            config['subcommands']['bedtools']['lmod']
        ))
        load_modules(config['subcommands']['bedtools']['lmod'])
    guide2flank_cmd = '{cmd} flank -i stdin -g "{chrom_sizes}" -b {flanking}'.format(
        cmd=' '.join(config['subcommands']['bedtools']['cmd']),
        chrom_sizes=reference_files['chrom_sizes'],
        flanking=config['parameters']['flanking'],
    )
    flank2seq_cmd = '{cmd} getfasta -fi "{genome_fasta}" -bed stdin -fo stdout -nameOnly -tab'.format(
        cmd=' '.join(config['subcommands']['bedtools']['cmd']),
        genome_fasta=reference_files['fasta'],
    )
    guide_flank_cmd = f'set -o pipefail; {guide2flank_cmd} | {flank2seq_cmd}'
    _logger.debug(f'Guide flank subcommand: {guide_flank_cmd}')
    return Subcommand(guide_flank_cmd, shell=True)

def bed_intersect_cmd(file_b, config):
    """Return a subcommand that reads in a bed file (a) and outputs the
    intersection with a second bed file (b).
    """
    if config['subcommands']['bedtools']['lmod']:
        _logger.debug('LMOD modules for guide flank subcommand: {}'.format(
            config['subcommands']['bedtools']['lmod']
        ))
        load_modules(config['subcommands']['bedtools']['lmod'])
    bedtools_intersect_cmd = [*config['subcommands']['bedtools']['cmd'],
                              'intersect', '-a', 'stdin', '-b', file_b]
    _logger.debug(f'BED intersect subcommand: {" ".join(bedtools_intersect_cmd)}')
    return Subcommand(*bedtools_intersect_cmd)
