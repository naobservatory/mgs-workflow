/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { DUSTMASKER_FASTQ_GZIPPED } from "../../../modules/local/dustmasker"
include { MINIMAP2_ONT as MINIMAP2_HV } from "../../../modules/local/minimap2"
include { SAMTOOLS_KEEP_AS_SAM } from "../../../modules/local/samtools"
include { CONCAT_GROUP_SINGLE as CONCAT_GROUP } from "../../../modules/local/concatGroup"
include { MERGE_SAM } from "../../../modules/local/samtools"

/***********
| WORKFLOW |
***********/

workflow EXTRACT_ONT_VIRAL_READS {
    take:
        reads
        group_ch
        grouping
        minimap2_hv_index
    main:
        // Grouping files
        if (grouping){
            grouped_ch = group_ch.join(reads, by: 0)
            .map { sample, group, reads -> tuple(sample, reads, group) }
            .groupTuple(by: 2)
            // Split into multi-sample and single-sample groups
            multi_sample_groups = grouped_ch.filter { it[0].size() > 1 }
            single_sample_groups = grouped_ch.filter { it[0].size() == 1 }
                .map { samples, read_list, group -> tuple(group, [read_list[0]]) }
            grouped_ch = CONCAT_GROUP(multi_sample_groups).mix(single_sample_groups)
        } else {
            grouped_ch = reads
        }
        // Drop non-complex reads
        masked_ch = DUSTMASKER_FASTQ_GZIPPED(grouped_ch)

        // Identify HV reads
        minimap2_hv_sam_ch = MINIMAP2_HV(masked_ch, minimap2_hv_index, "hv")
        minimap2_hv_sam_ch = SAMTOOLS_KEEP_AS_SAM(minimap2_hv_sam_ch, "hv")

        // Merging doesn't yet work, to fix.
        merged_sam_ch = MERGE_SAM(minimap2_hv_sam_ch.sam.collect(), "hv")

    emit:
        sam = merged_sam_ch.merged_sam
}

