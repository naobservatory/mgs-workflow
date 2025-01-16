/**************************************************
| SUBWORKFLOW: CONCATENATE & SUBSET RAW MGS READS |
**************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/
include { COUNT_READS } from "../../../modules/local/countReads"
include { MERGE_TSVS } from "../../../modules/local/mergeTsvs"

/***********
| WORKFLOW |
***********/

workflow COUNT_TOTAL_READS {
    take:
        samplesheet_ch
    main:
        read_counts_ch = COUNT_READS(samplesheet_ch)
        all_read_counts_ch = read_counts_ch.collect()
        read_counts_file = MERGE_TSVS(all_read_counts_ch, "read_counts")
    emit:
        read_counts = read_counts_file
}
