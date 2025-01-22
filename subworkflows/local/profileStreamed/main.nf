/**********************************************
| SUBWORKFLOW: HIGH-LEVEL TAXONOMIC PROFILING |
**********************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BBDUK_STREAMED as BBDUK } from "../../../modules/local/bbduk"
include { TAXONOMY_STREAMED as TAXONOMY_RIBO } from "../../../subworkflows/local/taxonomyStreamed"
include { TAXONOMY_STREAMED as TAXONOMY_NORIBO } from "../../../subworkflows/local/taxonomyStreamed"
include { ADD_FIXED_COLUMN as ADD_KRAKEN_RIBO } from "../../../modules/local/addFixedColumn"
include { ADD_FIXED_COLUMN as ADD_BRACKEN_RIBO } from "../../../modules/local/addFixedColumn"
include { ADD_FIXED_COLUMN as ADD_KRAKEN_NORIBO } from "../../../modules/local/addFixedColumn"
include { ADD_FIXED_COLUMN as ADD_BRACKEN_NORIBO } from "../../../modules/local/addFixedColumn"
include { CONCATENATE_TSVS as CONCATENATE_KRAKEN } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS as CONCATENATE_BRACKEN } from "../../../modules/local/concatenateTsvs"

include { MERGE_TAXONOMY_RIBO } from "../../../modules/local/mergeTaxonomyRibo"

/****************
| MAIN WORKFLOW |
****************/

workflow PROFILE_STREAMED {
    take:
        reads_ch
        kraken_db_ch
        ref_dir
        min_kmer_fraction
        k
        bbduk_suffix
        bracken_threshold
        single_end
    main:
        // Separate ribosomal reads
        ribo_path = "${ref_dir}/results/ribo-ref-concat.fasta.gz"
        ribo_ch = BBDUK(reads_ch, ribo_path, min_kmer_fraction, k, bbduk_suffix, !single_end)
        // Run taxonomic profiling separately on ribo and non-ribo reads
        tax_ribo_ch = TAXONOMY_RIBO(ribo_ch.match, kraken_db_ch, "D", bracken_threshold, single_end)
        tax_noribo_ch = TAXONOMY_NORIBO(ribo_ch.nomatch, kraken_db_ch, "D", bracken_threshold, single_end)
        // Add ribosomal status to output TSVs
        kr_ribo = ADD_KRAKEN_RIBO(tax_ribo_ch.kraken_reports, "ribosomal", "TRUE", "kraken_reports_ribo")
        kr_noribo = ADD_KRAKEN_NORIBO(tax_ribo_ch.kraken_reports, "ribosomal", "FALSE", "kraken_reports_noribo")
        br_ribo = ADD_BRACKEN_RIBO(tax_ribo_ch.bracken, "ribosomal", "TRUE", "bracken_ribo")
        br_noribo = ADD_BRACKEN_NORIBO(tax_ribo_ch.bracken, "ribosomal", "FALSE", "bracken_noribo")
        // Concatenate output TSVs
        kr_out = CONCATENATE_KRAKEN(kr_ribo.output.combine(kr_noribo.output), "kraken_reports_merged")
        br_out = CONCATENATE_BRACKEN(br_ribo.output.combine(br_noribo.output), "bracken_reports_merged")
    emit:
        bracken = br_out
        kraken = kr_out
}
