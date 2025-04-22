// Add a column to a TSV with sample ID
process RUN_LCA {
    label "python"
    label "single"
    input:
        path(tax_node_dmp)
        path(blast_tsv)
        val(group_column)
        val(taxid_column)
        val(output)
    output:
        path("${output}"), emit: output
        path("input_${blast_tsv}"), emit: input
    shell:
        '''
        run_lca.py !{tax_node_dmp} !{blast_tsv} !{group_column} !{taxid_column} !{output}
        # Link input files for testing
        ln -s !{blast_tsv} input_!{blast_tsv}
        '''
}
