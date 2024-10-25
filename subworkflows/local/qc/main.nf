/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { FASTQC } from "../../../modules/local/fastqc" addParams(cpus: "${params.fastqc_cpus}", mem: "${params.fastqc_mem}")
include { FALCO } from "../../../modules/local/falco" addParams(cpus: "${params.fastqc_cpus}", mem: "${params.fastqc_mem}")
include { MULTIQC } from "../../../modules/local/multiqc"
include { SUMMARIZE_MULTIQC_SINGLE } from "../../../modules/local/summarizeMultiqcSingle"

/***********
| WORKFLOW |
***********/

workflow QC {
    take:
        reads
        stage_label
    main:
        //fastqc_ch = FASTQC(reads)
        falco_ch = FALCO(reads)
        multiqc_ch = MULTIQC(stage_label, falco_ch.data.collect().ifEmpty([]))
        process_ch = SUMMARIZE_MULTIQC_SINGLE(multiqc_ch.data)
    emit:
        qc = process_ch
}
