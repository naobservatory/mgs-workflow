/**********************************************
| SUBWORKFLOW: HIGH-LEVEL TAXONOMIC PROFILING |
**********************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SUBSET_READS_PAIRED_TARGET; SUBSET_READS_PAIRED_TARGET as SUBSET_READS_PAIRED_TARGET_GROUP } from "../../../modules/local/subsetReads"
include { BBDUK } from "../../../modules/local/bbduk"
include { TAXONOMY as TAXONOMY_RIBO } from "../../../subworkflows/local/taxonomy"
include { TAXONOMY as TAXONOMY_NORIBO } from "../../../subworkflows/local/taxonomy"
include { MERGE_TAXONOMY_RIBO } from "../../../modules/local/mergeTaxonomyRibo"
include { CONCAT_GROUP } from "../../../modules/local/concatGroup"

/****************
| MAIN WORKFLOW |
****************/

workflow PROFILE {
    take:
        reads_ch
        group_ch
        kraken_db_ch
        n_reads
        ref_dir
        min_kmer_fraction
        k
        bbduk_suffix
        grouping
    main:
        // Randomly subset reads to target number
        subset_ch = SUBSET_READS_PAIRED_TARGET(reads_ch, n_reads, "fastq")
        if (grouping){
            // Join samplesheet with trimmed_reads and update fastq files
            subset_group_ch = group_ch.join(subset_ch, by: 0)
            .map { sample, group, reads -> tuple(sample, reads[0], reads[1], group) }
            .groupTuple(by: 3)
            // Split into multi-sample groups, these need to be subsetted to target number
            multi_sample_groups = subset_group_ch.filter { it[0].size() > 1 }
            // These are already subsetted to target number
            single_sample_groups = subset_group_ch.filter { it[0].size() == 1 }
                .map { samples, fwd_list, rev_list, group -> tuple(group, [fwd_list[0], rev_list[0]]) }
            // Concatenate multi-sample groups
            grouped_samples = CONCAT_GROUP(multi_sample_groups)
            // Randomly subset multi-sample groups to target number
            subset_grouped_ch = SUBSET_READS_PAIRED_TARGET_GROUP(grouped_samples, n_reads, "fastq")
            // Mix with subsetted multi-sample group with already subsetted single-sample groups
            grouped_ch = subset_grouped_ch.mix(single_sample_groups)
        } else {
            grouped_ch = subset_ch
        }
        // Separate ribosomal reads
        ribo_path = "${ref_dir}/results/ribo-ref-concat.fasta.gz"
        ribo_ch = BBDUK(grouped_ch, ribo_path, min_kmer_fraction, k, bbduk_suffix)
        // Run taxonomic profiling separately on ribo and non-ribo reads
        tax_ribo_ch = TAXONOMY_RIBO(ribo_ch.fail, kraken_db_ch, false, "D")
        tax_noribo_ch = TAXONOMY_NORIBO(ribo_ch.reads, kraken_db_ch, false, "D")
        // Merge ribo and non-ribo outputs
        kr_ribo = tax_ribo_ch.kraken_reports.collectFile(name: "kraken_reports_ribo.tsv.gz")
        kr_noribo = tax_noribo_ch.kraken_reports.collectFile(name: "kraken_reports_noribo.tsv.gz")
        br_ribo = tax_ribo_ch.bracken.collectFile(name: "bracken_reports_ribo.tsv.gz")
        br_noribo = tax_noribo_ch.bracken.collectFile(name: "bracken_reports_noribo.tsv.gz")
        merge_ch = MERGE_TAXONOMY_RIBO(kr_ribo, kr_noribo, br_ribo, br_noribo)
    emit:
        bracken = merge_ch.bracken
        kraken = merge_ch.kraken
}
