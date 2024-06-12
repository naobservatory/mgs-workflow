process MULTIQC {
    label "single"
    label "MultiQC"
    input:
        val(stage_label)
        path("*")
    output:
        path("multiqc_report.html"), emit: report
        tuple val(stage_label), path("multiqc_data"), emit: data
    shell:
        '''
        multiqc .
        '''
}
