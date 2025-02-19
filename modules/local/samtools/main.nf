process MERGE_SAM {
    label "samtools"

    input:
        path(sam_files)
        val prefix
    output:
        path("${prefix}_alignments.sam"), emit: merged_sam

    shell:
        '''
        # Merge SAM files and automatically add RG tags based on filenames
        samtools merge -r -O sam "!{prefix}_alignments.sam" !{sam_files}
        '''
}