/***********************************************************
| SUBWORKFLOW: ADAPTER & QUALITY TRIMMING OF RAW MGS READS |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { QC } from "../../../subworkflows/local/qc"
if (params.ont) {
    include { FILTLONG as FILTER_READS } from "../../../modules/local/filtlong"
} else {
    include { FASTP as FILTER_READS } from "../../../modules/local/fastp"
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
        single_end
    main:
        if (params.ont) {
            filter_ch = FILTER_READS(reads_ch)
        } else {
            filter_ch = FILTER_READS(reads_ch, adapter_path, single_end)
        }
        qc_ch = QC(filter_ch.reads, fastqc_cpus, fastqc_mem, stage_label, single_end)
    emit:
        reads = filter_ch.reads
        qc = qc_ch.qc
}