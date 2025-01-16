process COUNT_READS {
    label "coreutils_gzip_gawk"
    label "single_cpu_32GB_memory"
    input:
        tuple val(sample), path(reads)
    output:
        path("${sample}_read_count.tsv"), emit: output
        tuple val(sample), path("${sample}_reads_in_{1,2}.fastq.gz"), emit: input
    shell:
        '''
        # Count reads in file
        COUNT=$(zcat !{reads[0]} | wc -l | awk '{print $1/4}')
        # Add header
        echo -e "sample\ttotal_read_count" > !{sample}_read_count.tsv
        # Add sample and count
        echo -e "!{sample}\t$COUNT" >> !{sample}_read_count.tsv
        # Link output to input for tests
        ln -s !{reads[0]} !{sample}_reads_in_1.fastq.gz
        ln -s !{reads[1]} !{sample}_reads_in_2.fastq.gz
        '''
}
