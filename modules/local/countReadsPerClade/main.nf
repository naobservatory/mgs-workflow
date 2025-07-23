// Count reads per clade using taxonomic LCA assignments
// Note that if any reads are assigned to taxids not in the taxonomy, these reads will not be counted
process COUNT_READS_PER_CLADE {
    label "python"
    label "single"
    input:
        tuple val(sample), path(reads_tsv) // read table with LCA assignments and duplicate marking
        path(taxdb) // taxonomy database tsv with taxid and parent_taxid
    output:
        tuple val(sample), path("${sample}_clade_counts.tsv"), emit: output // tsv of clade counts per taxid
        tuple val(sample), path("input_${reads_tsv}"), emit: input
    script:
        """
        count_reads_per_clade.py --reads ${reads_tsv} --taxdb ${taxdb} --output ${sample}_clade_counts.tsv
        # Link input files for testing
        ln -s ${reads_tsv} input_${reads_tsv}
        """
}

