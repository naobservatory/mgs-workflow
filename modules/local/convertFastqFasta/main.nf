// Convert a single FASTQ file (interleaved or single-end) into FASTA format
// TODO: Expand to work on non-gzipped files

process CONVERT_FASTQ_FASTA {
    label "single"
    label "seqtk"
    input:
        tuple val(sample), path(fastq)
    output:
        tuple val(sample), path("${sample}_converted.fasta.gz"), emit: output
        tuple val(sample), path("${sample}_in.fastq.gz"), emit: input
    shell:
        '''
        # Perform conversion
        zcat !{fastq} | seqtk seq -a | gzip -c > !{sample}_converted.fasta.gz
        # Link input to output for testing
        ln -s !{fastq} !{sample}_in.fastq.gz
        '''
}
