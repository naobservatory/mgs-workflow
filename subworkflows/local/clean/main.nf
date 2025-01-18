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
    if (params.single_end) {
        include { FASTP_SINGLE as FILTER_READS } from "../../../modules/local/fastp"
    } else {
        include { FASTP_PAIRED as FILTER_READS } from "../../../modules/local/fastp"
    }
}
if (params.human_read_filtering) {
    include { MINIMAP2_ONT as MINIMAP2_HUMAN } from "../../../modules/local/minimap2"
    include { SAMTOOLS_FILTER } from "../../../modules/local/samtools"
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
        minimap2_human_index
    main:
        if (params.ont) {
            filter_ch = FILTER_READS(reads_ch)
            if (params.human_read_filtering) {
                minimap2_ch = MINIMAP2_HUMAN(filter_ch, minimap2_human_index, "human")
                filter_ch = SAMTOOLS_FILTER(minimap2_ch, "human")
            }
        } else {
            filter_ch = FILTER_READS(reads_ch, adapter_path)
        }
        qc_ch = QC(filter_ch.reads, fastqc_cpus, fastqc_mem, stage_label, single_end)
    emit:
        reads = filter_ch.reads
        qc = qc_ch.qc
}