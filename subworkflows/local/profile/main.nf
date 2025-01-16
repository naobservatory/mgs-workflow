/**********************************************
| SUBWORKFLOW: HIGH-LEVEL TAXONOMIC PROFILING |
**********************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/


if (params.single_end) {
    include { BBDUK_SINGLE as BBDUK } from "../../../modules/local/bbduk"
} else {
    include { BBDUK_PAIRED as BBDUK } from "../../../modules/local/bbduk"
}

include { BBDUK_HITS } from "../../../modules/local/bbduk"
include { TAXONOMY as TAXONOMY_RIBO } from "../../../subworkflows/local/taxonomy"
include { TAXONOMY as TAXONOMY_NORIBO } from "../../../subworkflows/local/taxonomy"
include { MERGE_TAXONOMY_RIBO } from "../../../modules/local/mergeTaxonomyRibo"

/****************
| MAIN WORKFLOW |
****************/

workflow PROFILE {
    take:
        reads_ch
        kraken_db_ch
        ref_dir
        min_kmer_fraction
        k
        bbduk_suffix
        single_end
    main:
        // Separate ribosomal reads
        ribo_path = "${ref_dir}/results/ribo-ref-concat.fasta.gz"
        ribo_ch = BBDUK(reads_ch, ribo_path, min_kmer_fraction, k, bbduk_suffix)
        // Run taxonomic profiling separately on ribo and non-ribo reads
        tax_ribo_ch = TAXONOMY_RIBO(ribo_ch.fail, kraken_db_ch, false, "D", single_end)
        tax_noribo_ch = TAXONOMY_NORIBO(ribo_ch.reads, kraken_db_ch, false, "D", single_end)
        // Merge ribo and non-ribo outputs
        kr_ribo = tax_ribo_ch.kraken_reports.collectFile(name: "kraken_reports_ribo.tsv.gz")
        kr_noribo = tax_noribo_ch.kraken_reports.collectFile(name: "kraken_reports_noribo.tsv.gz")
        br_ribo = tax_ribo_ch.bracken.collectFile(name: "bracken_reports_ribo.tsv.gz")
        br_noribo = tax_noribo_ch.bracken.collectFile(name: "bracken_reports_noribo.tsv.gz")
        merge_ch = MERGE_TAXONOMY_RIBO(kr_ribo, kr_noribo, br_ribo, br_noribo)
    emit:
        bracken = merge_ch.bracken
        kraken = merge_ch.kraken
        pre_qc = pre_qc_ch.qc
        post_qc = post_qc_ch.qc
}
