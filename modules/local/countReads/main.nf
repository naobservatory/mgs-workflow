process COUNT_READS {
    label "coreutils_gzip_gawk"
    label "single"
    input:
        tuple val(sample), path(reads)
        val(single_end)
    output:
        path("${sample}_read_count.tsv"), emit: output
        tuple val(sample), path("${sample}_reads_in.fastq.gz"), emit: input
    shell:
        '''
        # Count reads in file
        READS=!{single_end ? reads : reads[0]}
        COUNT=$(zcat ${READS} | wc -l | awk '{print $1/4}')
        COUNT_SINGLE=!{single_end ? '${COUNT}' : '$(expr ${COUNT} \\* 2)'}
        COUNT_PAIR=!{single_end ? 'NA' : '${COUNT}'}
        # Add header
        echo -e "sample\tn_reads_single\tn_read_pairs" > !{sample}_read_count.tsv
        # Add sample and count
        echo -e "!{sample}\t${COUNT_SINGLE}\t${COUNT_PAIR}" >> !{sample}_read_count.tsv
        # Link output to input for tests
        ln -s ${READS} !{sample}_reads_in.fastq.gz
        '''
}
