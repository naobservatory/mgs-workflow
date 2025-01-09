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
    shell:
        '''
        # Define input/output
        i=!{reads}
        o=!{sample}_!{suffix}_minimap2.sam
        ref=!{contaminant_ref}
        minimap2 -a ${ref} ${i} -o ${o}
        '''
} 