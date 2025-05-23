// Process tabular output from VSEARCH clustering into a unified output format
process PROCESS_VSEARCH_CLUSTER_OUTPUT {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(db)
        val n_clusters // Return representative IDs of the N largest clusters
        val output_prefix // Column name prefix for output DB
    output:
        tuple val(sample), path("${sample}_vsearch_tab.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_vsearch_ids.txt"), emit: ids
        tuple val(sample), path("input_${db}"), emit: input
    script:
        def prefix_string = output_prefix == "" ? "" : "-p ${output_prefix}"
        """
        in=${db}
        out_db=${sample}_vsearch_tab.tsv.gz
        out_id=${sample}_vsearch_ids.txt
        par="-n ${n_clusters} ${prefix_string}"
        process_vsearch_cluster_output.py \${par} \${in} \${out_db} \${out_id}
        ln -s \${in} input_\${in}
        """
}
