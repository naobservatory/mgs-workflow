process COUNT_DUPLICATE_READS {
    label "single"
    container "ubuntu:latest"
    
    input:
        tuple val(sample), path(tsv)
    
    output:
        tuple val(sample), path("${sample}_duplicate_reads.tsv")

    script:
    """
    tsv_processor "${tsv}" "${sample}_duplicate_reads.tsv" 0
    """
}