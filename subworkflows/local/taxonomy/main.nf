/***********************************************************
| SUBWORKFLOW: TAXONOMIC PROFILING WITH KRAKEN AND BRACKEN |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBMERGE } from "../../../modules/local/bbmerge"
include { CLUMPIFY_SINGLE } from "../../../modules/local/clumpify"

/***********
| WORKFLOW |
***********/

workflow TAXONOMY { // todo: rename
    take:
        reads_ch
    main:
        // Prepare reads
        joined_ch = BBMERGE(reads_ch)
        // Deduplicate reads
        dedup_ch = CLUMPIFY_SINGLE(joined_ch)
    emit:
        joined_reads = dedup_ch.reads
}
