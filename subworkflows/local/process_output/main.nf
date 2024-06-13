/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MERGE_TSVS } from "../../../modules/local/mergeTsvs"
include { SUMMARIZE_COMPOSITION } from "../../../modules/local/summarizeComposition"

/***********
| WORKFLOW |
***********/

workflow PROCESS_OUTPUT {
    take:
        multiqc_raw
        multiqc_cleaned
        multiqc_dedup
        multiqc_ribo_initial
        multiqc_ribo_secondary
        bracken_full
        bracken_pre
        bracken_post
        classify_dedup_subset
    main:
        // Collate MultiQC outputs
        multiqc_basic_ch = Channel.of(multiqc_raw[0], multiqc_cleaned[0], multiqc_dedup[0], multiqc_ribo_initial[0], multiqc_ribo_secondary[0])
        multiqc_adapter_ch = Channel.of(multiqc_raw[1], multiqc_cleaned[1], multiqc_dedup[1], multiqc_ribo_initial[1], multiqc_ribo_secondary[1])
        multiqc_qual_base_ch = Channel.of(multiqc_raw[2], multiqc_cleaned[2], multiqc_dedup[2], multiqc_ribo_initial[2], multiqc_ribo_secondary[2])
        multiqc_qual_seq_ch = Channel.of(multiqc_raw[3], multiqc_cleaned[3], multiqc_dedup[3], multiqc_ribo_initial[3], multiqc_ribo_secondary[3])
        // Merge MultiQC outputs
        basic_out_ch = MERGE_TSVS(multiqc_basic_ch.collect().ifEmpty([]))
        adapter_out_ch = MERGE_TSVS(multiqc_adapter_ch.collect().ifEmpty([]))
        qual_base_out_ch = MERGE_TSVS(multiqc_qual_base_ch.collect().ifEmpty([]))
        qual_seq_out_ch = MERGE_TSVS(multiqc_qual_seq_ch.collect().ifEmpty([]))
        // Summarize taxonomic composition
        tax_comp_full = SUMMARIZE_COMPOSITION(bracken_full, basic_out_ch, "ribo_secondary", 1)
        tax_comp_pre = SUMMARIZE_COMPOSITION(bracken_pre, basic_out_ch, "cleaned", classify_dedup_subset)
        tax_comp_pre = SUMMARIZE_COMPOSITION(bracken_post, basic_out_ch, "dedup", classify_dedup_subset)
    emit:
        basic = basic_out_ch
        adapt = adapter_out_ch
        qbase = qual_base_out_ch
        qseqs = qual_seq_out_ch
        composition_full = tax_comp_full
        composition_pre = tax_comp_pre
        composition_post = tax_comp_post
}
