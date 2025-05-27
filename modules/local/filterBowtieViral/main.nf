// modules/local/filterBowtieViral.nf
process FILTER_BOWTIE_VIRAL {
    label "samtools_coreutils"
    label "single"
    input:
        tuple val(sample), path(original_sam), path(reads_fastq)
    output:
        tuple val(sample), path("${sample}_noncontaminant.sam"), emit: sam
    script:
        """
        zcat ${reads_fastq} | grep "^@" -  | cut -d " " -f1 | sed 's/^@//' > temp_id.txt
        samtools view -h -N temp_id.txt -o ${sample}_noncontaminant.sam ${original_sam}
        """
}
