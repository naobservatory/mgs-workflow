// Copy a file to a new location with a custom path
process FASTQ_TO_FASTA {
    label "base"
    label "single"
    input:
        tuple val(sample), path(gzipped_reads)
    output:
        path("${sample}_R{1,2}.fasta.gz"), emit: fastas
    script:
    """
    # Convert R1
    zcat ${gzipped_reads[0]} | \
    awk 'NR%4==1{printf ">%s\\n", substr(\$0,2)}
         NR%4==2{print}' | \
    gzip > ${sample}_R1.fasta.gz
    
    # Convert R1
    zcat ${gzipped_reads[1]} | \
    awk 'NR%4==1{printf ">%s\\n", substr(\$0,2)}
         NR%4==2{print}' | \
    gzip > ${sample}_R2.fasta.gz

    """
}
