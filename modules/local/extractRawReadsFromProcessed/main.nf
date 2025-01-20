// Given paired & subset processed reads, extract the equivalent raw reads
process EXTRACT_RAW_READS_FROM_PROCESSED {
    label "seqtk"
    label "single"
    input:
        tuple val(sample), path(reads_processed), path(reads_raw)
        val(name)
    output:
        tuple val(sample), path("${sample}_${name}_{1,2}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_${name}_in_{processed,raw_1,raw_2}.fastq.gz"), emit: input
    shell:
        '''
        # Extract read IDs from filtered files
        seqtk comp !{reads_processed[0]} | cut -f 1 > read_ids.txt
        # Extract matching reads from raw files
        out_prefix="!{sample}_!{name}"
        seqtk subseq !{reads_raw[0]} read_ids.txt | gzip -c > ${out_prefix}_1.fastq.gz
        seqtk subseq !{reads_raw[1]} read_ids.txt | gzip -c > ${out_prefix}_2.fastq.gz
        # Link input to output for testing
        ln -s !{reads_processed} !{sample}_!{name}_in_processed.fastq.gz
        ln -s !{reads_raw[0]} !{sample}_!{name}_in_raw_1.fastq.gz
        ln -s !{reads_raw[1]} !{sample}_!{name}_in_raw_2.fastq.gz
        '''
}
