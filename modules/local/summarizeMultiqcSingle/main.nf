// Extract MultiQC data into more usable forms
process SUMMARIZE_MULTIQC_SINGLE {
    label "R"
    label "single"
    input:
        tuple val(stage), path(multiqc_data)
    output:
        tuple path("*_qc_basic_stats.tsv.gz"), path("*_qc_adapter_stats.tsv.gz"), path("*_qc_quality_base_stats.tsv.gz"), path("*_qc_quality_sequence_stats.tsv.gz")
    shell:
        '''
        summarize-multiqc-single.R -i !{multiqc_data} -s !{stage} -o ${PWD}
        '''
}
