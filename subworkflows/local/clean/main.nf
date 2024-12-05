/***********************************************************
| SUBWORKFLOW: ADAPTER & QUALITY TRIMMING OF RAW MGS READS |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc"
include { FASTP } from "../../../modules/local/fastp"

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
        single_end
    main:
        fastp_ch = FASTP(reads_ch, adapter_path, single_end)
        qc_ch = QC(fastp_ch.reads, fastqc_cpus, fastqc_mem, stage_label, single_end)
    emit:
        reads = fastp_ch.reads
        qc = qc_ch.qc
}