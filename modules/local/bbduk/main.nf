// Streamed version (interleaved or single-end input and output)
process BBDUK {
    label "small"
    label "BBTools"
    input:
        tuple val(sample), path(reads) // Interleaved or single-end
        path(contaminant_ref)
        val(params_map) // min_kmer_fraction, k, suffix, interleaved
    output:
        tuple val(sample), path("${sample}_${params_map.suffix}_bbduk_nomatch.fastq.gz"), emit: nomatch
        tuple val(sample), path("${sample}_${params_map.suffix}_bbduk_match.fastq.gz"), emit: match
        tuple val(sample), path("${sample}_${params_map.suffix}_bbduk.stats.txt"), emit: log
        tuple val(sample), path("${sample}_${params_map.suffix}_in.fastq.gz"), emit: input
    shell:
        '''
        suffix="!{params_map.suffix}"
        # Define input/output
        op=!{sample}_${suffix}_bbduk_nomatch.fastq.gz
        of=!{sample}_${suffix}_bbduk_match.fastq.gz
        stats=!{sample}_${suffix}_bbduk.stats.txt
        ref=!{contaminant_ref}
        il=!{params_map.interleaved ? 't' : 'f'}
        io="in=stdin.fastq ref=${ref} out=${op} outm=${of} stats=${stats} interleaved=${il}"
        # Define parameters
        par="minkmerfraction=!{params_map.min_kmer_fraction} k=!{params_map.k} t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        zcat !{reads} | bbduk.sh ${io} ${par}
        # Link input to output for testing
        ln -s !{reads} !{sample}_${suffix}_in.fastq.gz
        '''
}

// Streamed version of BBDUK_HITS that returns an interleaved file
// Uses minkmerhits instead of minkmerfraction
process BBDUK_HITS_INTERLEAVE {
    label "small"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        path(contaminant_ref)
        val(params_map) // min_kmer_hits, k, suffix
    output:
        tuple val(sample), path("${sample}_${params_map.suffix}_bbduk_in_{1,2}.fastq.gz"), emit: input
        tuple val(sample), path("${sample}_${params_map.suffix}_bbduk_pass.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_${params_map.suffix}_bbduk_fail.fastq.gz"), emit: fail
        tuple val(sample), path("${sample}_${params_map.suffix}_bbduk.stats.txt"), emit: log
    shell:
        '''
        suffix="!{params_map.suffix}"
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        op=!{sample}_${suffix}_bbduk_pass.fastq.gz
        of=!{sample}_${suffix}_bbduk_fail.fastq.gz
        stats=!{sample}_${suffix}_bbduk.stats.txt
        ref=!{contaminant_ref}
        io="in=stdin.fastq ref=${ref} out=${op} outm=${of} stats=${stats}"
        # Define parameters
        par="minkmerhits=!{params_map.min_kmer_hits} k=!{params_map.k} interleaved=t t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        # Execute
        paste <(zcat !{reads[0]} | paste - - - - ) <(zcat !{reads[1]} | paste - - - -) | tr "\t" "\n" | bbduk.sh ${io} ${par}
        # Move inputs for testing
        ln -s ${in1} !{sample}_${suffix}_bbduk_in_1.fastq.gz
        ln -s ${in2} !{sample}_${suffix}_bbduk_in_2.fastq.gz
        '''
}
