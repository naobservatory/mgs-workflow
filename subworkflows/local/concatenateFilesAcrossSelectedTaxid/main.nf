/*
Given a channel of 2-tuples, each containing a combined group/species label and
a file path, carry out the following operations:
- Replace the group/species label in each tuple with a group label only
- Group all files sharing a group label
- Concatenate all files within each group into a single file
The result is a restructured channel of 2-tuples, each containing a group label and
a concatenated file path.
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { CONCATENATE_FILES_BY_EXTENSION } from "../../../modules/local/concatenateFilesByExtension"

/***********
| WORKFLOW |
***********/

workflow CONCATENATE_FILES_ACROSS_SELECTED_TAXID {
    take:
        files // Channel of 2-tuples, each containing a combined group/species label and a file path
        filename_suffix // Suffix for output filenames
    main:
        // 1. Restructure channel to replace group/species label with group label only
        split_label_ch = files.map{
            label, path ->
                def pattern = /^(.*?)_(\d+)$/
                def matcher = (label =~ pattern)
                if (!matcher) {
                    def msg = "Group label doesn't match required pattern: ${label}, ${path}, ${pattern}"
                    throw new IllegalArgumentException(msg)
                }
                [matcher[0][1], path]
        }
        // 2. Group elements by group label and concatenate
        regrouped_label_ch = split_label_ch.groupTuple()
        regrouped_concat_ch = CONCATENATE_FILES_BY_EXTENSION(regrouped_label_ch, filename_suffix).output
    emit:
        output = regrouped_concat_ch
        test_input = files
}
