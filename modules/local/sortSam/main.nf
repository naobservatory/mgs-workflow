process SORT_SAM {
    label "pysam"
    input:
        tuple val(sample), path(sam)
    output:
        tuple val(sample), path("${sample}_sorted.sam"), emit: sam
        tuple val(sample), path("${sample}_sorted_nosq.sam"), emit: sam_nosq
    shell:
        '''
        smart_sort_sam.py !{sam} !{sample}_sorted.sam
        samtools view -h -G 12 !{sample}_sorted.sam | grep -v '^@SQ' - > !{sample}_sorted_nosq.sam
        '''
}
