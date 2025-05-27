process SORT_SAM {
    label "python"
    input:
        tuple val(sample), path(sam)
    output:
        tuple val(sample), path("${sample}_sorted.sam"), emit: sam
    shell:
        '''
        smart_sort_sam.py !{sam} !{sample}_sorted.sam
        '''
}
