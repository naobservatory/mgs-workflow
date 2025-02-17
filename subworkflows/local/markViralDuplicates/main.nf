/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SORT_TSV as SORT_INPUT } from "../../../modules/local/sortTsv"
include { SORT_TSV as SORT_GROUPS } from "../../../modules/local/sortTsv"
include { JOIN_TSVS } from "../../../modules/local/joinTsvs"
include { PARTITION_TSV } from "../../../modules/local/partitionTsv"
include { CONCATENATE_TSVS_LABELED } from "../../../modules/local/concatenateTsvs"

/***********
| WORKFLOW |
***********/

workflow MARK_VIRAL_DUPLICATES {
    take:
        groups // Labeled viral hit TSVs partitioned by group
    main:
        // 1. Sort inputs by sample ID
        input_ch = input_files.map{ label, input, groups -> tuple(label, input) }
        groups_ch = input_files.map{ label, input, groups -> tuple(label, groups) }
        input_sorted_ch = SORT_INPUT(input_ch, "sample").sorted
        groups_sorted_ch = SORT_GROUPS(groups_ch, "sample").sorted
        combined_sorted_ch = input_sorted_ch.combine(groups_sorted_ch, by: 0)
        // 2. Add group information to hit TSVs
        joined_ch = JOIN_TSVS(combined_sorted_ch, "sample", "strict", "input").output
        // 3. Partition each TSV by group ID
        partitioned_ch = PARTITION_TSV(joined_ch, "group").output
        // 4. Restructure channel so all files with the same group ID are together
        partitioned_ch.view()
        // First rearrange each element from [sample, [paths]] to [[group1, path1], [group2, path2], ...]
        partitioned_tuple_ch = partitioned_ch.map{
            sample, filepaths ->
                def pathlist = (filepaths instanceof List) ? filepaths : [filepaths]
                pathlist.collect { path ->
                    println(path)
                    def matcher = (path.last() =~ /^partition_(.*?)_${sample}_input_strict_joined_sample\.tsv\.gz$/)
                    if (!matcher) {
                        def msg = "Filename doesn't match required pattern: ${sample}, ${path}, ${path.last()}"
                        throw new IllegalArgumentException(msg)
                    }
                    [matcher[0][1], path]
                }
            }
        // Then rearrange channel to [[group1, [paths]], [group2, [paths]], ...]
        partitioned_grouped_ch = partitioned_tuple_ch.map{it[0]}.groupTuple()
        // 5. Concatenate TSVs for each group
        concat_ch = CONCATENATE_TSVS_LABELED(partitioned_grouped_ch, "grouped").output
    emit:
        groups = concat_ch
        test_in = input_files
        test_join = joined_ch
        test_part = partitioned_ch
        test_grps = partitioned_grouped_ch
}
