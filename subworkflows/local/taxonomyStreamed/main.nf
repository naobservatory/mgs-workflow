/***********************************************************
| SUBWORKFLOW: TAXONOMIC PROFILING WITH KRAKEN AND BRACKEN |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBMERGE_STREAMED as BBMERGE } from "../../../modules/local/bbmerge"
include { JOIN_FASTQ_STREAMED as JOIN_FASTQ } from "../../../modules/local/joinFastq"
include { KRAKEN_STREAMED as KRAKEN } from "../../../modules/local/kraken"
include { LABEL_KRAKEN_REPORTS_STREAMED as LABEL_KRAKEN_REPORTS } from "../../../modules/local/labelKrakenReports"
include { CONCATENATE_TSVS_STREAMED as CONCATENATE_KRAKEN_REPORTS } from "../../../modules/local/concatenateTsvs"
include { BRACKEN2 as BRACKEN } from "../../../modules/local/bracken"
include { LABEL_BRACKEN_REPORTS_STREAMED as LABEL_BRACKEN_REPORTS } from "../../../modules/local/labelBrackenReports"
include { CONCATENATE_TSVS_STREAMED as CONCATENATE_BRACKEN_REPORTS } from "../../../modules/local/concatenateTsvs"
include { SUMMARIZE_BBMERGE_STREAMED as SUMMARIZE_BBMERGE } from "../../../modules/local/summarizeBBMerge"

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
            single_read_ch = JOIN_FASTQ(merged_ch.reads, false).reads
            // Summarize the merged elements
            summarize_bbmerge_ch = SUMMARIZE_BBMERGE(merged_ch.reads).summary
        }
        // Run Kraken and munge reports
        kraken_ch = KRAKEN(single_read_ch, kraken_db_ch)
        kraken_label_ch = LABEL_KRAKEN_REPORTS(kraken_ch.report)
        kraken_merge_ch = CONCATENATE_KRAKEN_REPORTS(kraken_label_ch.report.collect().ifEmpty([]), "kraken_reports")
        // Run Bracken and munge reports
        bracken_ch = BRACKEN(kraken_ch.report, kraken_db_ch, classification_level) // NB: Not streamed
        bracken_label_ch = LABEL_BRACKEN_REPORTS(bracken_ch)
        bracken_merge_ch = CONCATENATE_BRACKEN_REPORTS(bracken_label_ch.report.collect().ifEmpty([]), "bracken_reports")
    emit:
        input_reads = reads_ch
        single_reads = single_read_ch
        bbmerge_summary = summarize_bbmerge_ch
        kraken_output = kraken_ch.output
        kraken_reports = kraken_merge_ch.output
        bracken = bracken_merge_ch.output
}
