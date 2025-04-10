// Subset a FASTA or FASTQ file to specific IDs using seqtk subseq
process SUBSEQ_FASTN {
    label "seqkit"
    label "single"
    input:
        tuple val(sample), path(fastn), path(ids)
    output:
        tuple val(sample), path("subseq_${fastn}"), emit: output
        tuple val(sample), path("input_${fastn}"), path("input_${ids}"),  emit: input
    shell:
        '''
        seqkit grep -f !{ids} !{fastn} | !{fastn.toString().endsWith(".gz") ? 'gzip -c' : 'cat'} > subseq_!{fastn}
        ln -s !{fastn} input_!{fastn}
        ln -s !{ids} input_!{ids}
        '''
}
