@Grab(group='org.codehaus.groovy', module='groovy-yaml', version='3.0.16')
import groovy.yaml.YamlBuilder

process PRIMER_DESIGN_CONFIG {
    label 'process_low'
    publishDir enabled: false

    input:
        val fasta
        val chrom_sizes
        val coding_bed

    output:
        path "pd_config.yaml", emit: pd_config

    script:
         def config = [
            references: [
                genome: [
                    fasta: "${fasta}",
                    coding_bed: [
                        CDS: "${coding_bed}"
                    ],
                    chrom_sizes: "${chrom_sizes}"
                ]
            ],
            ngs: [
                adapters: params._ngs
            ]
        ]
    
        YamlBuilder yaml = new YamlBuilder()
        yaml config
        def config_yaml = yaml.toString()

        """
        echo '${config_yaml}' > pd_config.yaml
        """
}
