// Process tabular output from VSEARCH clustering into a unified output format
process PROCESS_VSEARCH_CLUSTER_OUTPUT {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(db)
        val n_clusters // Return representative IDs of the N largest cluster
    output:
        tuple val(sample), path("${sample}_vsearch_tab.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_vsearch_ids.txt"), emit: ids
        tuple val(sample), path("input_${db}"), emit: input
    shell:
        '''
        in=!{db}
        out_db=!{sample}_vsearch_tab.tsv.gz
        out_id=!{sample}_vsearch_ids.txt
        process_vsearch_cluster_output.py -n !{n_clusters} ${in} ${out_db} ${out_id}
        ln -s ${in} input_${in}
        '''
}
