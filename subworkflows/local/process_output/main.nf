/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MERGE_TSVS as MERGE_MULTIQC_BASIC } from "../../../modules/local/mergeTsvs" addParams(name: "qc_basic_stats")
include { MERGE_TSVS as MERGE_MULTIQC_ADAPT } from "../../../modules/local/mergeTsvs" addParams(name: "qc_adapter_stats")
include { MERGE_TSVS as MERGE_MULTIQC_QBASE } from "../../../modules/local/mergeTsvs" addParams(name: "qc_quality_base_stats")
include { MERGE_TSVS as MERGE_MULTIQC_QSEQS } from "../../../modules/local/mergeTsvs" addParams(name: "qc_quality_sequence_stats")
include { SUMMARIZE_COMPOSITION as SUMMARIZE_COMPOSITION_FULL } from "../../../modules/local/summarizeComposition"
include { SUMMARIZE_COMPOSITION as SUMMARIZE_COMPOSITION_PRE } from "../../../modules/local/summarizeComposition"
include { SUMMARIZE_COMPOSITION as SUMMARIZE_COMPOSITION_POST } from "../../../modules/local/summarizeComposition"

/***********
| WORKFLOW |
***********/

workflow PROCESS_OUTPUT {
    take:
        multiqc_ch
        bracken_full
        bracken_pre
        bracken_post
        classify_dedup_subset
    main:
        // Collate MultiQC outputs
        multiqc_basic_ch = multiqc_ch.map{ it[0] }.collect().ifEmpty([])
        multiqc_adapt_ch = multiqc_ch.map{ it[1] }.collect().ifEmpty([])
        multiqc_qbase_ch = multiqc_ch.map{ it[2] }.collect().ifEmpty([])
        multiqc_qseqs_ch = multiqc_ch.map{ it[3] }.collect().ifEmpty([])
        // Merge MultiQC outputs
        basic_out_ch = MERGE_MULTIQC_BASIC(multiqc_basic_ch)
        adapt_out_ch = MERGE_MULTIQC_ADAPT(multiqc_adapt_ch)
        qbase_out_ch = MERGE_MULTIQC_QBASE(multiqc_qbase_ch)
        qseqs_out_ch = MERGE_MULTIQC_QSEQS(multiqc_qseqs_ch)
        // Summarize taxonomic composition
        tax_comp_full = SUMMARIZE_COMPOSITION_FULL(bracken_full, basic_out_ch, "ribo_secondary", 1)
        tax_comp_pre = SUMMARIZE_COMPOSITION_PRE(bracken_pre, basic_out_ch, "cleaned", classify_dedup_subset)
        tax_comp_post = SUMMARIZE_COMPOSITION_POST(bracken_post, basic_out_ch, "dedup", classify_dedup_subset)
    emit:
        basic = basic_out_ch
        adapt = adapt_out_ch
        qbase = qbase_out_ch
        qseqs = qseqs_out_ch
        composition_full = tax_comp_full
        composition_pre = tax_comp_pre
        composition_post = tax_comp_post
}
