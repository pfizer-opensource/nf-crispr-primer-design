//
// Subworkflow to deploy a reference genome for the crispr primer design tool
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDOPS_GTF2BED     } from '../../modules/nf-core/bedops/gtf2bed/main'
include { PREP_REFERENCE     } from '../../modules/local/prep-reference'

include { getGenomeAttribute } from './utils_nfcore_nf-crispr-primer-design_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO DEPLOY A REFERENCE GENOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEPLOY_REFERENCE {

    take:
    genome // string: Reference genome ID

    main:
    // Directory with reference files for this genome
    genome_dir = file("${params.genomes_dir}/${genome}")

    // Map of genome reference source files
    ref_sources = [
        fasta: file(getGenomeAttribute('fasta')),
        faidx: file(getGenomeAttribute('fasta') + '.fai'),
        bowtie2_1: file(getGenomeAttribute('bowtie2') + '/genome.1.bt2'),
        bowtie2_2: file(getGenomeAttribute('bowtie2') + '/genome.2.bt2'),
        bowtie2_3: file(getGenomeAttribute('bowtie2') + '/genome.3.bt2'),
        bowtie2_4: file(getGenomeAttribute('bowtie2') + '/genome.4.bt2'),
        bowtie2_r1: file(getGenomeAttribute('bowtie2') + '/genome.rev.1.bt2'),
        bowtie2_r2: file(getGenomeAttribute('bowtie2') + '/genome.rev.2.bt2'),
        gtf: file(getGenomeAttribute('gtf'))
    ]

    // Map of genome reference files on local storage
    ref_files = ref_sources.collectEntries { key, source ->
        [ key, genome_dir / "${source.name}" ]
    }
    // Additional files that are generated from other reference files
    ref_files += [
        chrom_sizes: genome_dir / "genome.chrom.sizes",
        genes_bed: genome_dir / "${ref_files.gtf.baseName}.bed",
        coding_bed: genome_dir / "coding.bed"
    ]

    // Check if all reference files are already in local storage
    if (ref_files.every { _, f -> f.exists() }) {
        // Collect bowtie2 reference files into a single map entry
        ref_files.bowtie2 = ref_files.findAll {k, _ -> k.startsWith("bowtie2")}.collect {it.value}
    }
    else {
        // Reference files are not present locally

        // Generate a genes BED file from the genes GTF file
        BEDOPS_GTF2BED([[:], ref_sources.gtf])
        BEDOPS_GTF2BED.out.bed
            .map { _meta, genes_bed -> genes_bed }
            .set { genes_bed }

        // Get the reference files, process them, and store them locally
        prep_ref_out = PREP_REFERENCE(
            genome,
            ref_sources.fasta,
            ref_sources.faidx,
            ref_sources.findAll {k, _ -> k.startsWith("bowtie2")}.collect {it.value},
            ref_sources.gtf,
            genes_bed
        )

        // Use the genome directory and reference files from the module output
        genome_dir = prep_ref_out.genome_dir
        ref_files = prep_ref_out
    }

    emit:
    genome_dir  = genome_dir
    fasta       = ref_files.fasta
    faidx       = ref_files.faidx
    chrom_sizes = ref_files.chrom_sizes
    bowtie2     = ref_files.bowtie2
    gtf         = ref_files.gtf
    genes_bed   = ref_files.genes_bed
    coding_bed  = ref_files.coding_bed
}
