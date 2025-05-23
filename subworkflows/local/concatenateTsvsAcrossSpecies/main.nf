/*
Given a channel of 2-tuples, each containing a combined group/species label and
a TSV path, carry out the following operations:
- Add a column to each TSV containing the group/species label
- Replace the group/species label in each tuple with a group label only
- Concatenate all TSVs sharing a group label into a single TSV
- Add a column to each new TSV containing the group label
The result is a restructured channel of 2-tuples, each containing a group label and
a TSV path, with the latter file now containing two additional columns: one for the
group/species label and one for the group label.
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SELECT_TSV_COLUMNS as DROP_GROUP_COLUMN } from "../../../modules/local/selectTsvColumns"
include { SELECT_TSV_COLUMNS as DROP_GROUP_SPECIES_COLUMN } from "../../../modules/local/selectTsvColumns"
include { ADD_SAMPLE_COLUMN as LABEL_GROUP_SPECIES } from "../../../modules/local/addSampleColumn"
include { CONCATENATE_TSVS_LABELED } from "../../../modules/local/concatenateTsvs"
include { ADD_SAMPLE_COLUMN as LABEL_GROUP } from "../../../modules/local/addSampleColumn"

/***********
| WORKFLOW |
***********/

workflow CONCATENATE_TSVS_ACROSS_SPECIES {
    take:
        tsvs // Channel of 2-tuples, each containing a combined group/species label and a TSV path
        filename_suffix // Suffix for output filenames
    main:
        // 0. Remove existing group and group_species columns, if present
        drop_group_ch = DROP_GROUP_COLUMN(tsvs, "group", "drop").output
        drop_group_species_ch = DROP_GROUP_SPECIES_COLUMN(drop_group_ch, "group_species", "drop").output
        // 1. Add a group/species label column to each TSV
        label_ch = LABEL_GROUP_SPECIES(drop_group_species_ch, "group_species", "group_species").output
        // 2. Restructure channel to replace group/species label with group label only
        split_label_ch = label_ch.map{
            label, path ->
                def pattern = /^(.*?)_(\d+)$/
                def matcher = (label =~ pattern)
                if (!matcher) {
                    def msg = "Group label doesn't match required pattern: ${label}, ${path}, ${pattern}"
                    throw new IllegalArgumentException(msg)
                }
                [matcher[0][1], path]
        }
        // 3. Group elements by group label and concatenate
        regrouped_label_ch = split_label_ch.groupTuple()
        regrouped_concat_ch = CONCATENATE_TSVS_LABELED(regrouped_label_ch, filename_suffix).output
        // 4. Add a group label column to each new TSV
        regrouped_ch = LABEL_GROUP(regrouped_concat_ch, "group", "group").output
    emit:
        output = regrouped_ch
        test_input = tsvs
}
