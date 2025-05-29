// Consolidated viral SAM filtering: removes contaminant reads, applies pair-based score threshold filtering, adds missing mates for UP reads, and sorts output
process FILTER_VIRAL_SAM {
    label "python"
    label "single"
    input:
        tuple val(sample), path(sam), path(contaminant_reads)
        val(score_threshold)
    output:
        tuple val(sample), path("${sample}_viral_filtered.sam"), emit: sam
    script:
        """
        # Extract contaminant read IDs from FASTQ
        zcat ${contaminant_reads} | grep "^@" | cut -d " " -f1 | sed 's/^@//' > contaminant_ids.txt
        
        # Run the consolidated viral SAM filtering
        filter_viral_sam.py ${sam} contaminant_ids.txt ${sample}_viral_filtered.sam ${score_threshold}
        """
}