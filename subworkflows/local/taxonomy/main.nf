/***********************************************************
| SUBWORKFLOW: TAXONOMIC PROFILING WITH KRAKEN AND BRACKEN |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_PAIRED } from "../../../modules/local/subsetReads" addParams(suffix: "fastq")
include { BBMERGE } from "../../../modules/local/bbmerge"
include { JOIN_FASTQ } from "../../../modules/local/joinFastq"
include { CLUMPIFY_SINGLE } from "../../../modules/local/clumpify"
include { KRAKEN } from "../../../modules/local/kraken" addParams(mem: "${params.kraken_memory}")
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
            subset_ch = SUBSET_READS_PAIRED(reads_ch, params.read_fraction)
        }
        // Prepare reads
        merged_ch = BBMERGE(subset_ch)
        joined_ch = JOIN_FASTQ(merged_ch.reads)
        // Deduplicate reads (if applicable)
        if ( params.dedup_rc ){
            dedup_ch = CLUMPIFY_SINGLE(joined_ch)
        } else {
            dedup_ch = joined_ch
        }
        // Run Kraken and munge reports
        kraken_ch = KRAKEN(joined_ch, kraken_db_ch)
        kraken_label_ch = LABEL_KRAKEN_REPORTS(kraken_ch.report)
        kraken_merge_ch = MERGE_KRAKEN_REPORTS(kraken_label_ch.collect().ifEmpty([]))
        // Run Bracken and munge reports
        bracken_ch = BRACKEN(kraken_ch.report, kraken_db_ch, params.classification_level)
        bracken_label_ch = LABEL_BRACKEN_REPORTS(bracken_ch)
        bracken_merge_ch = MERGE_BRACKEN(bracken_label_ch.collect().ifEmpty([]))
    emit:
        kraken_output = kraken_ch.output
        kraken_reports = kraken_merge_ch
        bracken = bracken_merge_ch
}
