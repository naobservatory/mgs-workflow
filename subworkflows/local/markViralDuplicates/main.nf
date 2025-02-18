/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { MARK_ALIGNMENT_DUPLICATES } from "../../../modules/local/markAlignmentDuplicates"

/***********
| WORKFLOW |
***********/

workflow MARK_VIRAL_DUPLICATES {
    take:
        groups // Labeled viral hit TSVs partitioned by group
        deviation // Maximum alignment deviation that qualifies as a duplicate
    main:
        dup_ch = MARK_ALIGNMENT_DUPLICATES(groups, deviation).output
    emit:
        dup = dup_ch
        test_in = groups
}
