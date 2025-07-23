// Count reads per clade using taxonomic LCA assignments
// Note that if any reads are assigned to taxids not in the taxonomy, these reads will not be counted
process COUNT_READS_PER_CLADE {
    label "python"
    label "single"
    input:
        tuple val(sample), path(reads_tsv) // read table with columns: seq_id, prim_align_dup_exemplar, aligner_taxid_lca
        path(taxdb) // taxonomy database tsv with columns: taxid, parent_taxid
    output:
        tuple val(sample), path("${sample}_clade_counts.tsv.gz"), emit: output // gzipped tsv with columns: taxid, parent_taxid, reads_direct_total, reads_direct_dedup, reads_clade_total, reads_clade_dedup
        tuple val(sample), path("input_${reads_tsv}"), emit: input
    script:
        """
        count_reads_per_clade.py --reads ${reads_tsv} --taxdb ${taxdb} --output ${sample}_clade_counts.tsv.gz
        # Link input files for testing
        ln -s ${reads_tsv} input_${reads_tsv}
        """
}

