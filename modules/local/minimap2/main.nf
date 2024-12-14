// Detection and removal of contaminant reads
process MINIMAP2 {
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
        minimap2 -t $(nproc) -ax map-ont ${ref} ${i} -o ${o}
        '''
}

# Reference to use: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz