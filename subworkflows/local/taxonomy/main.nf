/***********************************************************
| SUBWORKFLOW: TAXONOMIC PROFILING WITH KRAKEN AND BRACKEN |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

if (params.read_type == "paired_end") {
    include { SUBSET_READS_PAIRED as SUBSET_READS } from "../../../modules/local/subsetReads" addParams(suffix: "fastq")
} else if (params.read_type == "single_end") {
    include { SUBSET_READS_SINGLE as SUBSET_READS } from "../../../modules/local/subsetReads" addParams(suffix: "fastq")
}
include { BBMERGE } from "../../../modules/local/bbmerge" // probalby  skippable in single read version
include { SUMMARIZE_BBMERGE } from "../../../modules/local/summarizeBBMerge" // probalby  skippable in single read version
include { SUMMARIZE_DEDUP } from "../../../modules/local/summarizeDedup" // probalby  no change needed in single read version
include { CLUMPIFY_PAIRED } from "../../../modules/local/clumpify" // already has a single read version.
include { JOIN_FASTQ } from "../../../modules/local/joinFastq" // probably not needed in single read version
include { CLUMPIFY_SINGLE } from "../../../modules/local/clumpify" // already has a single read version.
include { KRAKEN } from "../../../modules/local/kraken" addParams(mem: "${params.kraken_memory}") // probs not changes needed
include { LABEL_KRAKEN_REPORTS } from "../../../modules/local/labelKrakenReports"
include { MERGE_TSVS as MERGE_KRAKEN_REPORTS } from "../../../modules/local/mergeTsvs" addParams(name: "kraken_reports")
include { MERGE_TSVS as MERGE_BRACKEN } from "../../../modules/local/mergeTsvs" addParams(name: "bracken_reports")
include { BRACKEN } from "../../../modules/local/bracken"
include { LABEL_BRACKEN_REPORTS } from "../../../modules/local/labelBrackenReports"

/***********
| WORKFLOW |
***********/

workflow TAXONOMY {
    take:
        reads_ch
        kraken_db_ch
    main:
        // Subset reads (if applicable)
        if ( params.read_fraction == 1 ){
            subset_ch = reads_ch
        } else {
            subset_ch = SUBSET_READS(reads_ch, params.read_fraction)
        }

        if (params.read_type == "paired_end") {
            // Deduplicate reads (if applicable)
            if ( params.dedup_rc ){
                paired_dedup_ch = CLUMPIFY_PAIRED(subset_ch)
            } else {
                paired_dedup_ch = subset_ch
            }
            // Prepare reads
            merged_ch = BBMERGE(paired_dedup_ch)
            // Only want to summarize the merged elements
            summarize_bbmerge_ch = SUMMARIZE_BBMERGE(merged_ch.reads.map{sample, files -> [sample, files[0]]})
            joined_ch = JOIN_FASTQ(merged_ch.reads)
            // Deduplicate reads (if applicable)
            if ( params.dedup_rc ){
                dedup_ch = CLUMPIFY_SINGLE(joined_ch)
            } else {
                dedup_ch = joined_ch
            }

        } else if (params.read_type == "single_end") {
            // Deduplicate reads (if applicable)

            if (params.dedup_rc) {
                dedup_ch = CLUMPIFY_SINGLE(subset_ch)
            } else {
                dedup_ch = subset_ch
            }
        }

        // Summarize last of the output
        summarize_dedup_ch = SUMMARIZE_DEDUP(dedup_ch)

        // Run Kraken and munge reports
        kraken_ch = KRAKEN(dedup_ch, kraken_db_ch)
        kraken_label_ch = LABEL_KRAKEN_REPORTS(kraken_ch.report)
        kraken_merge_ch = MERGE_KRAKEN_REPORTS(kraken_label_ch.collect().ifEmpty([]))
        // Run Bracken and munge reports
        bracken_ch = BRACKEN(kraken_ch.report, kraken_db_ch, params.classification_level)
        bracken_label_ch = LABEL_BRACKEN_REPORTS(bracken_ch)
        bracken_merge_ch = MERGE_BRACKEN(bracken_label_ch.collect().ifEmpty([]))

        if (params.read_type == "paired_end") {
            emit:
                kraken_output = kraken_ch.output
                kraken_reports = kraken_merge_ch
                bracken = bracken_merge_ch
                bbmerge_summary = summarize_bbmerge_ch
                dedup_summary = summarize_dedup_ch
        }

        else if (params.read_type == "single_end") {
            emit:
                kraken_output = kraken_ch.output
                kraken_reports = kraken_merge_ch
                bracken = bracken_merge_ch
                dedup_summary = summarize_dedup_ch
        }
}
