// Process tabular output from VSEARCH clustering into a unified output format
process PROCESS_VSEARCH_CLUSTER_OUTPUT {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(db)
    output:
        tuple val(sample), path("${sample}_vsearch_tab.tsv.gz"), emit: output
        tuple val(sample), path("input_${db}"), emit: input
    shell:
        '''
        in=!{db}
        out=!{sample}_vsearch_tab.tsv.gz
        process_vsearch_cluster_output.py ${in} ${out}
        ln -s ${in} input_${in}
        '''
}
