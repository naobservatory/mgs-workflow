/**************************************************
| SUBWORKFLOW: CONCATENATE & SUBSET RAW MGS READS |
**************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/
include { COUNT_READS } from "../../../modules/local/countReads"
include { CONCATENATE_TSVS } from "../../../modules/local/concatenateTsvs"

/***********
| WORKFLOW |
***********/

workflow COUNT_TOTAL_READS {
    take:
        samplesheet_ch
    main:
        read_counts_ch = COUNT_READS(samplesheet_ch)
        all_read_counts_ch = read_counts_ch.output.collect()
        read_counts_file = CONCATENATE_TSVS(all_read_counts_ch, "read_counts")
    emit:
        read_counts = read_counts_file
}
