// Unleave an interleaved FASTQ file and save paired output files
// TODO: Implement streamed version

process UNLEAVE_FASTQ {
    label "small"
    label "BBTools"
    input:
        tuple val(sample), path(reads_interleaved)
    output:
        tuple val(sample), path("${sample}_unleaved_{1,2}.fastq.gz"), emit: output
        tuple val(sample), path("${sample}_in.fastq.gz"), emit: input
    shell:
        '''
        par="t=!{task.cpus} -Xmx!{task.memory.toGiga()}g"
        io="in=!{reads_interleaved} out=!{sample}_unleaved_1.fastq.gz out2=!{sample}_unleaved_2.fastq.gz"
        reformat.sh ${io} ${par}
        ln -s !{reads_interleaved} !{sample}_in.fastq.gz # Link input to output for testing
        '''
}
