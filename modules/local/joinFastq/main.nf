// Join and concatenate partially-merged FASTQ files into a single read file
process JOIN_FASTQ {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_joined.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_joined_in.fastq.gz"), emit: input
    shell:
        '''
        # Prepare to join unmerged read pairs
        om=!{reads[0]}
        ou1=!{reads[1]}
        ou2=!{reads[2]}
        oj=!{sample}_bbmerge_unmerged_joined.fastq.gz
        # Join unmerged read pairs
        join_fastq.py ${ou1} ${ou2} ${oj}
        # Concatenate single output file
        oo=!{sample}_joined.fastq.gz
        cat ${om} ${oj} > ${oo}
        # Link input reads for testing
        in2=!{sample}_joined_in.fastq.gz
        ln -s !{reads} ${in2}
        '''
}

// Join and concatenate partially-merged interleaved FASTQ files into a single read file
// TODO: Consider replacing with a Rust script for speed
process JOIN_FASTQ_STREAMED {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(reads) // Merged single reads, then unmerged interleaved
        val(debug)
    output:
        tuple val(sample), path("${sample}_joined.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_joined_in_{merged,unmerged}.fastq.gz"), emit: input
    shell:
        '''
        # Prepare to join unmerged read pairs
        om=!{reads[0]}
        ou=!{reads[1]}
        oj=!{sample}_bbmerge_unmerged_joined.fastq.gz
        # Join unmerged read pairs
        join_fastq_interleaved.py ${ou} ${oj} !{ debug ? "--debug" : "" }
        # Concatenate single output file
        oo=!{sample}_joined.fastq.gz
        cat ${om} ${oj} > ${oo}
        # Link input reads for testing
        im=!{sample}_joined_in_merged.fastq.gz
        iu=!{sample}_joined_in_unmerged.fastq.gz
        ln -s ${om} ${im}
        ln -s ${ou} ${iu}
        '''
}
