/*
Given a channel of 2-tuples, each containing a combined group/species label and
a FASTA path, carry out the following operations:
- Replace the group/species label in each tuple with a group label only
- Group all FASTA files sharing a group label 
- Concatenate all FASTA files within each group into a single FASTA file
The result is a restructured channel of 2-tuples, each containing a group label and
a concatenated FASTA path.
*/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { CONCATENATE_FASTN_LABELED } from "../../../modules/local/concatenateFastn"

/***********
| WORKFLOW |
***********/

workflow CONCATENATE_FASTA_ACROSS_SPECIES {
    take:
        fastas // Channel of 2-tuples, each containing a combined group/species label and a FASTA path
        filename_suffix // Suffix for output filenames
    main:
        // 1. Restructure channel to replace group/species label with group label only
        split_label_ch = fastas.map{
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
        regrouped_concat_ch = CONCATENATE_FASTN_LABELED(regrouped_label_ch, filename_suffix).output
    emit:
        output = regrouped_concat_ch
        test_input = fastas
}