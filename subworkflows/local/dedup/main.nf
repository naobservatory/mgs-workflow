/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc" addParams(fastqc_cpus: params.fastqc_cpus, fastqc_mem: params.fastqc_mem)
include { CLUMPIFY_PAIRED } from "../../../modules/local/clumpify"

/***********
| WORKFLOW |
***********/

workflow DEDUP {
    take:
        reads_ch
    main:
        dedup_ch = CLUMPIFY_PAIRED(reads_ch)
        qc_ch = QC(dedup_ch.reads, params.stage_label)
    emit:
        reads = dedup_ch.reads
        qc = qc_ch.out.qc
}
