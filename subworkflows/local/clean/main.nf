/***********************************************************
| SUBWORKFLOW: ADAPTER & QUALITY TRIMMING OF RAW MGS READS |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem)
include { CUTADAPT } from "../../../modules/local/cutadapt"
include { FASTP } from "../../../modules/local/fastp"

/***********
| WORKFLOW |
***********/

workflow CLEAN {
    take:
        reads_ch
        adapter_path
    main:
        adapt_ch = CUTADAPT(reads_ch, adapter_path)
        fastp_ch = FASTP(adapt_ch.reads, adapter_path)
        qc_ch = QC(fastp_ch.reads, params.stage_label)
    emit:
        reads = fastp_ch.reads
        qc = qc_ch.qc
}
