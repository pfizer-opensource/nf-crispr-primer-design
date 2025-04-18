/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { DEPLOY_REFERENCE       } from '../subworkflows/local/deploy-reference'

include { PRIMER_DESIGN_CONFIG   } from '../modules/local/primer-design-config'
include { PRIMER_DESIGN          } from '../modules/local/primer-design'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CRISPR_PRIMER_DESIGN {

    take:
    ch_input // channel: input file

    main:

    ch_versions = Channel.empty()

    // Deploy reference genome
    DEPLOY_REFERENCE(params.genome)

    // Generate config file for primer design tool
    PRIMER_DESIGN_CONFIG(
        DEPLOY_REFERENCE.out.fasta,
        DEPLOY_REFERENCE.out.chrom_sizes,
        DEPLOY_REFERENCE.out.coding_bed
    )

    // Run primer design tool
    PRIMER_DESIGN(
        Channel.of([:]).merge(ch_input),
        "genome",
        PRIMER_DESIGN_CONFIG.out.pd_config
    )
    ch_versions = ch_versions.mix(PRIMER_DESIGN.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nf-crispr-primer-design_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
