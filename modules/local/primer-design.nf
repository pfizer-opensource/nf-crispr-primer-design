process PRIMER_DESIGN {
    label 'process_low'
    label 'image_crispr_primer_design'

    input:
        tuple val(meta), path(guide_table)
        val reference

    output:
        path "${prefix}_CRISPR_primers.txt",     emit: primers
        path "${prefix}_CRISPR_primers_NGS.txt", emit: NGS_primers, optional: true
        path "${prefix}_amplicons.yaml",         emit: amplicons, optional: true
        path "${prefix}_CRISPR_primers.bed",     emit: bed, optional: true
        path "versions.yml",                     emit: versions

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        CRISPR_primer_design \\
            -i "${guide_table}" \\
            -o "${prefix}_CRISPR_primers.txt" \\
            -r "${reference}" \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            crispr_primer_design: \$(CRISPR_primer_design --version)
        END_VERSIONS
        """
}