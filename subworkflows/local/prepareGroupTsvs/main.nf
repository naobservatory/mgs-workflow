/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SORT_TSV as SORT_INPUT } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_GROUPS } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_JOINED_GROUPS } from "../../../modules/local/sortTsv"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { PARTITION_TSV } from "../../../modules/local/partitionTsv"
include { CONCATENATE_TSVS_LABELED } from "../../../modules/local/concatenateTsvs"

/***********
| WORKFLOW |
***********/

workflow PREPARE_GROUP_TSVS {
    take:
        input_files
    main:
        // 1. Sort inputs by sample ID
        input_ch = input_files.map{ label, input, _groups -> tuple(label, input) }
        groups_ch = input_files.map{ label, _input, groups -> tuple(label, groups) }
        input_sorted_ch = SORT_INPUT(input_ch, "sample").sorted
        groups_sorted_ch = SORT_GROUPS(groups_ch, "sample").sorted
        combined_sorted_ch = input_sorted_ch.combine(groups_sorted_ch, by: 0)
        // 2. Add group information to hit TSVs
        joined_ch = JOIN_TSVS(combined_sorted_ch, "sample", "inner", "input").output
        joined_sorted_ch = SORT_JOINED_GROUPS(joined_ch, "group").sorted
        // 3. Partition each TSV by group ID
        partitioned_ch = PARTITION_TSV(joined_sorted_ch, "group").output
        // 4. Restructure channel so all files with the same group ID are together
        // First rearrange each element from [sample, [paths]] to [[group1, path1], [group2, path2], ...]
        partitioned_flattened_ch = partitioned_ch.flatMap{
            sample, filepaths ->
                def pathlist = (filepaths instanceof List) ? filepaths : [filepaths]
                pathlist.collect { path ->
                    def filename = path.last()
                    def pattern = "^partition_(.*?)_sorted_group_${sample}_input_inner_joined_sample\\.tsv\\.gz\$"
                    def matcher = (filename =~ pattern)
                    if (!matcher) {
                        def msg = "Filename doesn't match required pattern: ${sample}, ${path}, ${path.last()}"
                        throw new IllegalArgumentException(msg)
                    }
                    [matcher[0][1], path]
                }
            }
        // Then rearrange channel to [[group1, [paths]], [group2, [paths]], ...]
        partitioned_grouped_ch = partitioned_flattened_ch.groupTuple()
        // 5. Concatenate TSVs for each group
        concat_ch = CONCATENATE_TSVS_LABELED(partitioned_grouped_ch, "grouped").output
    emit:
        groups = concat_ch
        test_in = input_files
        test_join = joined_ch
        test_part = partitioned_ch
        test_grps = partitioned_grouped_ch
}
