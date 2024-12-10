// Detection and removal of contaminant reads
process MINIMAP2 {
    label "large"
    label "minimap2"
    input:
        tuple val(sample), path(reads)
        path(contaminant_ref)
        // val(k)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_minimap2_contam.fastq.gz"), emit: sam
    shell:
        '''
        # Define input/output
        in=!{reads}
        op=!{sample}_!{suffix}_minimap2_pass.fastq.gz
        of=!{sample}_!{suffix}_minimap2_fail.fastq.gz

        ref=!{contaminant_ref}
        io="ref=${ref} in=${in} out=${op} outm=${of} stats=${stats}"
        var="-ax map-ont"
        # Define parameters
        # // par="t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        minimap2 ${io} ${var}
        '''
}