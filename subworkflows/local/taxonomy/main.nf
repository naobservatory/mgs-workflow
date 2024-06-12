/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_PAIRED } from "../modules/local/subsetReads"
include { BBMERGE } from "../modules/local/bbmerge"
include { CLUMPIFY_SINGLE } from "../modules/local/clumpify"
include { KRAKEN } from "../modules/local/kraken"
include { LABEL_KRAKEN_REPORTS } from "../modules/local/labelKrakenReports"
include { MERGE_TSVS } from "../modules/local/mergeTsvs"
include { BRACKEN } from "../modules/local/bracken"
include { LABEL_BRACKEN_REPORTS } from "../modules/local/labelBrackenReports"

/***********
| WORKFLOW |
***********/

workflow TAXONOMY {
    take:
        reads_ch
        kraken_db_ch
        read_fraction
        dedup_rc
        classification_level
    main:
        // Subset reads (if applicable)
        if ( read_fraction == 1 ){
            subset_ch = reads_ch
        } else {
            subset_ch = SUBSET_READS_PAIRED(reads_ch, read_fraction)
        }
        // Prepare reads
        merged_ch = BBMERGE(subset_ch)
        joined_ch = JOIN_FASTQ(merged_ch)
        // Deduplicate reads (if applicable)
        if ( dedup_rc ){
            dedup_ch = CLUMPIFY_SINGLE(joined_ch)
        } else {
            dedup_ch = joined_ch
        }
        // Run Kraken and munge reports
        kraken_ch = KRAKEN(joined_ch, kraken_db_ch)
        kraken_label_ch = LABEL_KRAKEN_REPORTS(kraken_ch.report)
        kraken_merge_ch = MERGE_TSVS(kraken_label_ch.collect().ifEmpty([]))
        // Run Bracken and munge reports
        bracken_ch = BRACKEN(kraken_ch.report, kraken_db_ch, classification_level)
        bracken_label_ch = LABEL_BRACKEN_REPORTS(bracken_ch)
        bracken_merge_ch = MERGE_TSVS(bracken_label_ch.collect().ifEmpty([]))
    emit:
        kraken_output = kraken_ch.output
        kraken_reports = kraken_merge_ch
        bracken = bracken_merge_ch
}
