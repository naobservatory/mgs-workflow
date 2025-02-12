// Detection and removal of contaminant reads, using indices created for ONT cDNA data
process MINIMAP2_ONT {
    label "large"
    label "minimap2"
    input:
        tuple val(sample), path(reads)
        path(contaminant_ref)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_minimap2.sam"), emit: sam
        tuple val(sample), path("${sample}_in.fastq.gz"), emit: input
    shell:
        '''
        # Define input
        o=!{sample}_!{suffix}_minimap2.sam
        ref=!{contaminant_ref}
        # Run minimap2
        zcat !{reads} | minimap2 -a ${ref} /dev/fd/0 > ${o}
        # Link input to output for testing
        ln -s !{reads} !{sample}_in.fastq.gz
        '''
}