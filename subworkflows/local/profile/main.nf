/**********************************************
| SUBWORKFLOW: HIGH-LEVEL TAXONOMIC PROFILING |
**********************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/


if (params.single_end) {
    include { SUBSET_READS_SINGLE_TARGET as SUBSET_READS_TARGET } from "../../../modules/local/subsetReads"
    include { BBDUK_SINGLE as BBDUK } from "../../../modules/local/bbduk"
    include { CONCAT_GROUP_SINGLE as CONCAT_GROUP } from "../../../modules/local/concatGroup"
    include { SUBSET_READS_SINGLE_TARGET; SUBSET_READS_SINGLE_TARGET as SUBSET_READS_TARGET_GROUP } from "../../../modules/local/subsetReads"
    include { FASTP_SINGLE as FASTP } from "../../../modules/local/fastp"
} else {
    include { SUBSET_READS_PAIRED_TARGET as SUBSET_READS_TARGET } from "../../../modules/local/subsetReads"
    include { SUBSET_READS_PAIRED_TARGET; SUBSET_READS_PAIRED_TARGET as SUBSET_READS_TARGET_GROUP } from "../../../modules/local/subsetReads"
    include { BBDUK_PAIRED as BBDUK } from "../../../modules/local/bbduk"
    include { CONCAT_GROUP_PAIRED as CONCAT_GROUP } from "../../../modules/local/concatGroup"
    include { FASTP_PAIRED as FASTP } from "../../../modules/local/fastp"
}

include { QC as PRE_ADAPTER_TRIM_QC } from "../../../subworkflows/local/qc"
include { QC as POST_ADAPTER_TRIM_QC } from "../../../subworkflows/local/qc"
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
        group_ch
        kraken_db_ch
        n_reads
        ref_dir
        min_kmer_fraction
        k
        bbduk_suffix
        grouping
        adapter_path
        fastqc_cpus
        fastqc_mem
        single_end
    main:


        // Randomly subset reads to target number
        subset_ch = SUBSET_READS_TARGET(reads_ch, n_reads, "fastq")

        if (grouping){
            // Join samplesheet with trimmed_reads and update fastq files
            if (single_end) {
                subset_group_ch = group_ch.join(subset_ch, by: 0)
                .map { sample, group, reads -> tuple(sample, reads, group) }
                .groupTuple(by: 2)
                // Single-sample groups are already subsetted to target number
                single_sample_groups = subset_group_ch.filter { it[0].size() == 1 }
                    .map { samples, read_list, group -> tuple(group, [read_list[0]]) }

            } else {
                subset_group_ch = group_ch.join(subset_ch, by: 0)
                .map { sample, group, reads -> tuple(sample, reads[0], reads[1], group) }
                .groupTuple(by: 3)
                single_sample_groups = subset_group_ch.filter { it[0].size() == 1 }
                    .map { samples, fwd_list, rev_list, group -> tuple(group, [fwd_list[0], rev_list[0]]) }
            }
            // Split into multi-sample groups, these need to be subsetted to target number
            multi_sample_groups = subset_group_ch.filter { it[0].size() > 1 }
            // Concatenate multi-sample groups
            grouped_samples = CONCAT_GROUP(multi_sample_groups)
            // Randomly subset multi-sample groups to target number
            subset_grouped_ch = SUBSET_READS_TARGET_GROUP(grouped_samples, n_reads, "fastq")
            // Mix with subsetted multi-sample group with already subsetted single-sample groups
            grouped_ch = subset_grouped_ch.mix(single_sample_groups)
        } else {
            grouped_ch = subset_ch
        }

        // Run FASTQC
        pre_qc_ch = PRE_ADAPTER_TRIM_QC(grouped_ch, fastqc_cpus, fastqc_mem, "pre_profile", single_end)

        // Call fastp adapter trimming
        fastp_ch = FASTP(grouped_ch, adapter_path)
        // Extract fastp trimmed reads
        trimmed_grouped_ch = fastp_ch.reads
        // Run FASTQC
        post_qc_ch = POST_ADAPTER_TRIM_QC(trimmed_grouped_ch, fastqc_cpus, fastqc_mem, "post_profile", single_end)

        // Separate ribosomal reads
        ribo_path = "${ref_dir}/results/ribo-ref-concat.fasta.gz"
        ribo_ch = BBDUK(trimmed_grouped_ch, ribo_path, min_kmer_fraction, k, bbduk_suffix)
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
