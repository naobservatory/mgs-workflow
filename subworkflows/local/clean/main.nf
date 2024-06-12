/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem)
include { CUTADAPT } from "../modules/local/cutadapt"
include { TRIMMOMATIC } from "../modules/local/trimmomatic"
include { FASTP } from "../modules/local/fastp"

/***********
| WORKFLOW |
***********/

workflow CLEAN {
    take:
        reads_ch
        adapter_path
    main:
        adapt_ch = CUTADAPT(reads_ch, adapter_path)
        trim_ch = TRIMMOMATIC(adapt_ch.reads, adapter_path)
        fastp_ch = FASTP(trim_ch.reads, adapter_path)
        qc_ch = QC(fastp_ch.reads, params.stage_label)
    emit:
        reads = fastp_ch.reads
        qc = qc_ch.out.qc
}
