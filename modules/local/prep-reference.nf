process PREP_REFERENCE {
    tag "${genome}"
    label 'process_low'
    publishDir enabled: false
    storeDir "${params.genomes_dir}"

    input:
        val genome
        path fasta
        path faidx
        path bowtie2
        path gtf
        path genes_bed

    output:
        path "${genome}",                    emit: genome_dir
        path "${genome}/${fasta}",           emit: fasta
        path "${genome}/${faidx}",           emit: faidx
        path "${genome}/genome.chrom.sizes", emit: chrom_sizes
        path "${genome}/*.bt2",              emit: bowtie2
        path "${genome}/${gtf}",             emit: gtf
        path "${genome}/${genes_bed}",       emit: genes_bed
        path "${genome}/coding.bed",         emit: coding_bed

    script:
        """
        # Per-genome reference directory
        mkdir -p "${genome}"

        # Copy reference files into genome reference directory, resolving symlinks
        cp -LR \\
            "${fasta}" \\
            "${faidx}" \\
            ${bowtie2} \\
            "${gtf}" \\
            "${genes_bed}" \\
            "${genome}/"

        # Generate chromosome sizes file (needed for bedtools) from faidx file
        cut -f1,2 "${faidx}" > "${genome}/genome.chrom.sizes"

        # Select the coding regions from the genes bed file
        awk '/protein_coding\tCDS/' < "${genes_bed}" > "${genome}/coding.bed"
        """
}