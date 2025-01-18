// Mask gzipped FASTA file with Dustmasker
process DUSTMASKER_FASTA_GZIPPED {
    label "BLAST"
    label "single"
    input:
        path(fasta_gzipped)
    output:
        path("masked.fasta.gz")
    shell:
        '''
        zcat -f !{fasta_gzipped} | dustmasker -out "masked.fasta" -outfmt fasta
        sed -i '/^>/!s/[a-z]/x/g' masked.fasta
        gzip masked.fasta
        '''
}

// Mask gzipped FASTQ file with Dustmasker
process DUSTMASKER_FASTQ_GZIPPED {
    label "BLAST"
    label "single"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_masked.fasta.gz")
    shell:
        '''
        # TURN FASTQ INTO FASTA
        zcat -f !{reads} | sed -n '1~4s/^@/>/p;2~4p' > !{sample}.fasta
        # MASK FASTA
        dustmasker -in !{sample}.fasta -out !{sample}_masked.fasta -outfmt fasta
        # Mark masked bases as N
        sed -i '/^>/!s/[a-z]/N/g' !{sample}_masked.fasta
        # GZIP FASTA
        gzip !{sample}_masked.fasta
        '''
}
