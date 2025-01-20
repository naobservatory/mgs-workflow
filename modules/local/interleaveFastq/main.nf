// Interleave paired FASTQ files into a single interleaved file
// TODO: Expand to work on non-gzipped files

process INTERLEAVE_FASTQ {
    label "single"
    label "coreutils_gzip_gawk"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_interleaved.fastq.gz"), emit: output
        tuple val(sample), path("${sample}_in_{1,2}.fastq.gz"), emit: input
    shell:
        '''
        # Perform interleaving
        paste <(zcat !{reads[0]} | paste - - - - ) <(zcat !{reads[1]} | paste - - - - ) | tr "\t" "\n" | gzip -c > "!{sample}_interleaved.fastq.gz"
        # Link input to output for testing
        ln -s !{reads[0]} !{sample}_in_1.fastq.gz
        ln -s !{reads[1]} !{sample}_in_2.fastq.gz
        '''
}
