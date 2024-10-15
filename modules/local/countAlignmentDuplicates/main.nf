process COUNT_ALIGNMENT_DUPLICATES {
    label "single"
    container "ubuntu:latest"
    
    input:
        tuple val(sample), path(tsv)
    
    output:
        path("${sample}_duplicate_reads.tsv")

    script:
    """
    find_duplicates "${tsv}" "${sample}_duplicate_reads.tsv" ${params.fuzzy_match}
    """
}