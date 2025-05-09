// Subset a FASTA or FASTQ file to specific IDs
process DOWNSAMPLE_FASTN_BY_ID {
    label "seqkit"
    label "single"
    input:
        tuple val(sample), path(fastn), path(ids)
    output:
        tuple val(sample), path("downsampled_${fastn}"), emit: output
        tuple val(sample), path("input_${fastn}"), path("input_${ids}"),  emit: input
    shell:
        '''
        seqkit grep -f !{ids} !{fastn} | seqkit rmdup | !{fastn.toString().endsWith(".gz") ? 'gzip -c' : 'cat'} > downsampled_!{fastn}
        ln -s !{fastn} input_!{fastn}
        ln -s !{ids} input_!{ids}
        '''
}
