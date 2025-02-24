// Pull out original clean reads from cleaned fastq, based on smaller FASTQ read ids
process PULLOUT_FASTQ {
    label "seqtk"
    label "single"
    input:
        tuple val(sample), path(mapped_fastq), path(cleaned_fastq)
    output:
        tuple val(sample), path("${sample}_cleaned_subset.fastq.gz"), emit: output
        tuple val(sample), path("${sample}_cleaned_in.fastq.gz"), path("${sample}_mapped_in.fastq.gz"), emit: input
    shell:
        '''
        set -euo pipefail

        # Get mapped read ids
        zcat !{mapped_fastq} | awk 'NR%4==1' | sed 's/^@//' > mapped_read_ids.txt

        # Create FASTQ file with cleaned, mapped reads
        seqtk subseq !{cleaned_fastq} mapped_read_ids.txt | gzip -c > "!{sample}_cleaned_subset.fastq.gz"

        # Link input to output for testing
        ln -s !{mapped_fastq} !{sample}_mapped_in.fastq.gz
        ln -s !{cleaned_fastq} !{sample}_cleaned_in.fastq.gz        '''
}
