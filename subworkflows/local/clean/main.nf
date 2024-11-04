/***********************************************************
| SUBWORKFLOW: ADAPTER & QUALITY TRIMMING OF RAW MGS READS |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem, stage_label: params.stage_label)
include { FASTP } from "../../../modules/local/fastp"

/***********
| WORKFLOW |
***********/

workflow CLEAN {
    take:
        reads_ch
        adapter_path
    main:
        fastp_ch = FASTP(reads_ch, adapter_path)
        qc_ch = QC(fastp_ch.reads)
    emit:
        reads = fastp_ch.reads
        qc = qc_ch.qc
}
