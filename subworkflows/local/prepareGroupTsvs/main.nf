/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

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
        // 1. Add group information to hit TSVs
        joined_ch = JOIN_TSVS(input_files, "sample", "strict", "input")
        // 2. Partition each TSV by group ID
        partitioned_ch = PARTITION_TSV(joined_ch, "group")
        // 3. Restructure channel so all files with the same group ID are together
        // First rearrange each element from [sample, [paths]] to [[group1, path1], [group2, path2], ...]
        partitioned_tuple_ch = partitioned_ch.output
            | map { sample, filepaths ->
                filepaths.map { path ->
                    def matcher = (path =~ /^${sample}_(.*?)_part\.tsv\.gz$/)
                    if (!matcher) {
                        def msg = "Filename doesn't match required pattern: ${sample}, ${path}"
                        throw new IllegalArgumentException(msg)
                    [matcher[0][1], fp]
                }
            }
        // Then rearrange channel to [[group1, [paths]], [group2, [paths]], ...]
        partitioned_grouped_ch = partitioned_tuple_ch.groupTuple()
        // 4. Concatenate TSVs for each group
        concat_ch = CONCATENATE_TSVS_LABELED(partitioned_grouped_ch, "grouped")
    emit:
        groups = concat_ch
        // TODO: Add more outputs for testing as needed
}
