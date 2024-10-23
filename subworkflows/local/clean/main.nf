/***********************************************************
| SUBWORKFLOW: ADAPTER & QUALITY TRIMMING OF RAW MGS READS |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem)
if (params.read_type == "single_end") {
    include { FASTP_SINGLE as FASTP } from "../../../modules/local/fastp"
} else if (params.read_type == "paired_end") {
    include { FASTP_PAIRED as FASTP } from "../../../modules/local/fastp"
}

/***********
| WORKFLOW |
***********/

workflow CLEAN {
    take:
        reads_ch
        adapter_path
    main:
        fastp_ch = FASTP(reads_ch, adapter_path)
        qc_ch = QC(fastp_ch.reads, params.stage_label)
    emit:
        reads = fastp_ch.reads
        qc = qc_ch.qc
}
