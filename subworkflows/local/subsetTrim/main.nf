/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

if (params.single_end) {
    include { SUBSET_READS_SINGLE_TARGET as SUBSET_READS_TARGET } from "../../../modules/local/subsetReads"
    include { CONCAT_GROUP_SINGLE as CONCAT_GROUP } from "../../../modules/local/concatGroup"
    include { SUBSET_READS_SINGLE_TARGET; SUBSET_READS_SINGLE_TARGET as SUBSET_READS_TARGET_GROUP } from "../../../modules/local/subsetReads"
    include { FASTP_SINGLE as FASTP } from "../../../modules/local/fastp"
} else {
    include { SUBSET_READS_PAIRED_TARGET as SUBSET_READS_TARGET } from "../../../modules/local/subsetReads"
    include { SUBSET_READS_PAIRED_TARGET; SUBSET_READS_PAIRED_TARGET as SUBSET_READS_TARGET_GROUP } from "../../../modules/local/subsetReads"
    include { CONCAT_GROUP_PAIRED as CONCAT_GROUP } from "../../../modules/local/concatGroup"
    include { FASTP_PAIRED as FASTP } from "../../../modules/local/fastp"
}

/***********
| WORKFLOW |
***********/

workflow SUBSET_TRIM {
    take:
      reads_ch
      group_ch
      n_reads
      grouping
      adapter_path
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

        // Call fastp adapter trimming
        fastp_ch = FASTP(grouped_ch, adapter_path)
    emit:
        subset_reads = grouped_ch
        trimmed_subset_reads = fastp_ch.reads
}
