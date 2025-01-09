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
        tuple val(sample), path("outfile_masked.fasta.gz")
    shell:
        '''
        # TODO: Make it so the files are named after the sample
        zcat -f !{reads} | sed -n '1~4s/^@/>/p;2~4p' > outfile.fasta
        dustmasker -in outfile.fasta -out outfile_masked.fasta -outfmt fasta
        sed -i '/^>/!s/[a-z]/x/g' outfile_masked.fasta
        gzip outfile_masked.fasta
        '''
}
