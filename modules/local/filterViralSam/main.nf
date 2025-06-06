// Consolidated viral SAM filtering: keep contaminant-free reads, applies pair-based score threshold filtering, adds missing mates for UP reads, and sorts output
process FILTER_VIRAL_SAM {
    label "pysam_biopython"
    label "single"
    input:
        tuple val(sample), path(sam), path(contaminant_free_reads)
        val(score_threshold)
    output:
        tuple val(sample), path("${sample}_viral_filtered.sam.gz"), emit: sam
    script:
        """
        sorted_fastq="${sample}_sorted_contaminant_free.fastq.gz"
        sorted_sam="${sample}_sorted_viral_filtered.sam.gz"
        outf="${sample}_viral_filtered.sam.gz"

        # Make sure fastq file is sorted
        zcat ${contaminant_free_reads} | paste - - - - | sort -k1,1 | tr '\\t' '\\n' | gzip -c > \${sorted_fastq}

        # Make sure SAM file is sorted
        zcat ${sam} | sort -t \$'\\t' -k1,1 | gzip -c > \${sorted_sam}

        # Run the consolidated viral SAM filtering
        filter_viral_sam.py \${sorted_sam} \${sorted_fastq} \${outf} ${score_threshold}
        """
}
