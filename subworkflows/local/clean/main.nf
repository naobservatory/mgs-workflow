/***********************************************************
| SUBWORKFLOW: ADAPTER & QUALITY TRIMMING OF RAW MGS READS |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc"
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
        fastqc_cpus
        fastqc_mem
        stage_label
        read_type
    main:
        fastp_ch = FASTP(reads_ch, adapter_path)
        qc_ch = QC(fastp_ch.reads, stage_label, read_type, fastqc_cpus, fastqc_mem)
    emit:
        reads = fastp_ch.reads
        qc = qc_ch.qc
}
