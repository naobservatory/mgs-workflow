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
        tuple val(sample), path("${sample}_masked.fasta.gz"), emit: masked
        tuple val(sample), path("${sample}_in.fastq.gz"), emit: input
    shell:
        '''
        set -e
        reads=!{reads}
        zcat -f ${reads} | \
            # Convert FASTQ to FASTA
            sed -n '1~4s/^@/>/p;2~4p' | \
            # Mask with Dustmasker
            dustmasker -in /dev/stdin -out /dev/stdout -outfmt fasta | \
            # Replace masked bases with N
            sed '/^>/!s/[a-z]/N/g' | \
            # Gzip FASTA
            gzip > !{sample}_masked.fasta.gz

        # Link input to output for testing
        ln -s ${reads} !{sample}_in.fastq.gz
        '''
}