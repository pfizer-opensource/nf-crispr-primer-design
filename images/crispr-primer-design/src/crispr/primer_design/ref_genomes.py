"""Reference genome handling."""
from pathlib import Path

from .exceptions import PrimerDesignGenomeError


def get_reference_files(references, genome, coding=None):
    if genome not in references:
        raise PrimerDesignGenomeError(
            'Invalid name of reference genome: {}'.format(genome)
        )
    reference_files = references[genome]
    try:
        fasta = Path(reference_files['fasta'])
    except KeyError:
        raise PrimerDesignGenomeError('Reference genome FASTA path required') from None
    if not fasta.exists():
        raise PrimerDesignGenomeError(
            'Reference genome FASTA file is not accessible: {}'.format(fasta)
        )
    if reference_files.get('genome_base'):
        genome_base = Path(reference_files['genome_base'])
    else:
        genome_base = fasta.parent / fasta.name.partition('.')[0]
        reference_files['genome_base'] = str(genome_base)
    if reference_files.get('chrom_sizes'):
        chrom_sizes = Path(reference_files['chrom_sizes'])
    else:
        chrom_sizes = genome_base.with_suffix('.chrom.sizes')
        reference_files['chrom_sizes'] = str(chrom_sizes)
    if not genome_base.with_suffix('.1.bt2').exists():
        raise PrimerDesignGenomeError(
            'Could not find bowtie2 index files at base: {}'.format(genome_base)
        )
    if not chrom_sizes.exists():
        raise PrimerDesignGenomeError(
            'Could not find chromesome sizes file: {}'.format(chrom_sizes)
        )
    if coding is not None:
        try:
            coding_bed = Path(reference_files['coding_bed'][coding])
        except KeyError:
            raise PrimerDesignGenomeError(
                'Invalid coding exon definition name: {}'.format(coding)
            )
        if not coding_bed.exists():
            raise PrimerDesignGenomeError(
                'Coding exon BED file is not accessible: {}'.format(coding_bed)
            )
    
    return reference_files

def get_ref_names(references):
    rv = {}
    for name, ref in sorted(references.items(), key=lambda k: k[1].priority or 0, reverse=True):
        rv[name] = sorted(ref.coding_bed, key=lambda x: (x == 'none', 1))
    return rv
