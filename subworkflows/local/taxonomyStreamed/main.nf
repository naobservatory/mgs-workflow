/***********************************************************
| SUBWORKFLOW: TAXONOMIC PROFILING WITH KRAKEN AND BRACKEN |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBMERGE_STREAMED as BBMERGE } from "../../../modules/local/bbmerge"
include { JOIN_FASTQ_STREAMED as JOIN_FASTQ } from "../../../modules/local/joinFastq"

include { SUMMARIZE_BBMERGE } from "../../../modules/local/summarizeBBMerge"
include { SUMMARIZE_DEDUP } from "../../../modules/local/summarizeDedup"
include { CLUMPIFY_PAIRED } from "../../../modules/local/clumpify"
include { CLUMPIFY_SINGLE } from "../../../modules/local/clumpify"
include { KRAKEN } from "../../../modules/local/kraken"
include { LABEL_KRAKEN_REPORTS } from "../../../modules/local/labelKrakenReports"
include { MERGE_TSVS as MERGE_KRAKEN_REPORTS } from "../../../modules/local/mergeTsvs"
include { MERGE_TSVS as MERGE_BRACKEN } from "../../../modules/local/mergeTsvs"
include { BRACKEN } from "../../../modules/local/bracken"
include { LABEL_BRACKEN_REPORTS } from "../../../modules/local/labelBrackenReports"

/***********
| WORKFLOW |
***********/

workflow TAXONOMY_STREAMED {
    take:
        reads_ch // Should be interleaved for paired-end data
        kraken_db_ch
        classification_level
        single_end
    main:
        if (single_end) {
            // No merging in single read version
            summarize_bbmerge_ch = Channel.empty()
            single_read_ch = reads_ch
        } else {
            // NB: No deduplication in streamed version
            // Merge and concatenate input reads
            merged_ch = BBMERGE(reads_ch)
            single_read_ch = JOIN_FASTQ(merged_ch.reads, false)
            // Only want to summarize the merged elements
            //summarize_bbmerge_ch = SUMMARIZE_BBMERGE(merged_ch.reads.map{sample, files -> [sample, files[0]]})
            //single_read_ch = JOIN_FASTQ(merged_ch.reads)
        }

//        // Deduplicate reads (if applicable)
//        if (dedup_rc) {
//                dedup_ch = CLUMPIFY_SINGLE(single_read_ch)
//            } else {
//                dedup_ch = single_read_ch
//        }
//
//        // Summarize last of the output
//        summarize_dedup_ch = SUMMARIZE_DEDUP(dedup_ch)
//
//        // Run Kraken and munge reports
//        kraken_ch = KRAKEN(dedup_ch, kraken_db_ch)
//        kraken_label_ch = LABEL_KRAKEN_REPORTS(kraken_ch.report)
//        kraken_merge_ch = MERGE_KRAKEN_REPORTS(kraken_label_ch.collect().ifEmpty([]), "kraken_reports")
//        // Run Bracken and munge reports
//        bracken_ch = BRACKEN(kraken_ch.report, kraken_db_ch, classification_level)
//        bracken_label_ch = LABEL_BRACKEN_REPORTS(bracken_ch)
//        bracken_merge_ch = MERGE_BRACKEN(bracken_label_ch.collect().ifEmpty([]), "bracken_reports")
    emit:
        test_out_reads = single_read_ch.reads
//        kraken_output = kraken_ch.output
//        kraken_reports = kraken_merge_ch
//        bracken = bracken_merge_ch
//        bbmerge_summary = summarize_bbmerge_ch
//        dedup_summary = summarize_dedup_ch
}
