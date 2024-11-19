// Extract paired MultiQC data into a more usable form
process SUMMARIZE_MULTIQC_PAIR {
    label "R"
    label "single"
    input:
        tuple val(stage), val(sample), path(multiqc_data)
        val(single_end)
    output:
        tuple path("${stage}_${sample}_qc_basic_stats.tsv.gz"), path("${stage}_${sample}_qc_adapter_stats.tsv.gz"), path("${stage}_${sample}_qc_quality_base_stats.tsv.gz"), path("${stage}_${sample}_qc_quality_sequence_stats.tsv.gz")
    shell:
        '''
        # Convert boolean to R format
        if [ "!{single_end}" = "true" ]; then
            single_end="TRUE"
        else
            single_end="FALSE"
        fi
        summarize-multiqc-pair.R -i !{multiqc_data} -s !{stage} -S !{sample} -r $single_end -o ${PWD}
        '''
}