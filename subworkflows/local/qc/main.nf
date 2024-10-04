/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

if (params.read_type == "single_end") {
    include { SUMMARIZE_MULTIQC_SINGLE as SUMMARIZE_MULTIQC } from "../../../modules/local/summarizeMultiqc"
} else if (params.read_type == "paired_end") {
    include { SUMMARIZE_MULTIQC_PAIRED as SUMMARIZE_MULTIQC } from "../../../modules/local/summarizeMultiqc"
}
include { FASTQC } from "../../../modules/local/fastqc" addParams(cpus: "${params.fastqc_cpus}", mem: "${params.fastqc_mem}")
include { MULTIQC } from "../../../modules/local/multiqc"

/***********
| WORKFLOW |
***********/

workflow QC {
    take:
        reads
        stage_label
    main:
        fastqc_ch = FASTQC(reads)
        multiqc_ch = MULTIQC(stage_label, fastqc_ch.zip.collect().ifEmpty([]))
        process_ch = SUMMARIZE_MULTIQC(multiqc_ch.data)
    emit:
        qc = process_ch
}