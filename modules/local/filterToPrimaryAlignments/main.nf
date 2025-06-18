// Filter the joined tsv of the LCA results with the Bowtie2 results with multiple 
// alignments to only keep the primary alignments for each read, such that we can 
// run the DOWNSTREAM workflow on the output of this pipeline.
process FILTER_TO_PRIMARY_ALIGNMENTS {
    label "python"
    label "single"
    input:
        tuple val(sample), path(joined_tsv)
    output:
        tuple val(sample), path("${sample}_primary_alignments.tsv.gz"), emit: output
        tuple val(sample), path("input_${joined_tsv}"), emit: input
    script:
        """
        # Set strict error handling
        set -o pipefail
        out=${sample}_primary_alignments.tsv.gz
        # Filter the joined tsv to only have the primary alignments
        zcat ${joined_tsv} | filter_to_primary_alignments.py -o \${out}
        # Link input to output for testing
        ln -s ${joined_tsv} input_${joined_tsv}
        """
}
