/***********************************************************
| SUBWORKFLOW: TAXONOMIC PROFILING WITH KRAKEN AND BRACKEN |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBMERGE } from "../../../modules/local/bbmerge"
include { JOIN_FASTQ } from "../../../modules/local/joinFastq"
include { SUMMARIZE_BBMERGE } from "../../../modules/local/summarizeBBMerge"
include { KRAKEN } from "../../../modules/local/kraken"
include { HEAD_TSV as HEAD_KRAKEN_REPORTS } from "../../../modules/local/headTsv"
include { ADD_SAMPLE_COLUMN as LABEL_KRAKEN_REPORTS } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS as CONCATENATE_KRAKEN_REPORTS } from "../../../modules/local/concatenateTsvs"
include { BRACKEN } from "../../../modules/local/bracken"
include { REHEAD_TSV as REHEAD_BRACKEN } from "../../../modules/local/reheadTsv"
include { ADD_SAMPLE_COLUMN as LABEL_BRACKEN } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS as CONCATENATE_BRACKEN } from "../../../modules/local/concatenateTsvs"

/***********
| WORKFLOW |
***********/

workflow TAXONOMY {
    take:
        reads_ch // Should be interleaved for paired-end data
        kraken_db_ch
        classification_level
        bracken_threshold
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
        kraken_headers = "pc_reads_total,n_reads_clade,n_reads_direct,n_minimizers_total,n_minimizers_distinct,rank,taxid,name"
        kraken_head_ch = HEAD_KRAKEN_REPORTS(kraken_ch.report, kraken_headers, "kraken_report")
        kraken_label_ch = LABEL_KRAKEN_REPORTS(kraken_head_ch.output, "sample", "kraken_report")
        kraken_combined_ch = kraken_label_ch.output.map{ sample, file -> file }.collect().ifEmpty([])
        kraken_merge_ch = CONCATENATE_KRAKEN_REPORTS(kraken_combined_ch, "kraken_reports")
        // Run Bracken and munge reports
        bracken_ch = BRACKEN(kraken_ch.report, kraken_db_ch, classification_level, bracken_threshold) // NB: Not streamed
        bracken_rehead_ch = REHEAD_BRACKEN(bracken_ch, "taxonomy_id,taxonomy_lvl", "taxid,rank")
        bracken_label_ch = LABEL_BRACKEN(bracken_rehead_ch.output, "sample", "bracken")
        bracken_combined_ch = bracken_label_ch.output.map{ sample, file -> file }.collect().ifEmpty([])
        bracken_merge_ch = CONCATENATE_BRACKEN(bracken_combined_ch, "bracken")
    emit:
        input_reads = reads_ch
        single_reads = single_read_ch
        bbmerge_summary = summarize_bbmerge_ch
        kraken_output = kraken_ch.output
        kraken_reports = kraken_merge_ch.output
        bracken = bracken_merge_ch.output
}
