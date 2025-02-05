// Concatenate multiple TSVs (streamed version with Python)
process CONCATENATE_TSVS {
    label "python"
    label "single"
    input:
        path(tsvs)
        val(name)
    output:
        path("${name}.tsv.gz"), emit: output
        path("${name}_in_0.tsv.gz"), emit: input
    shell:
        '''
        concatenate_tsvs.py -o !{name}.tsv.gz !{tsvs}
        ln -s !{tsvs[0]} !{name}_in_0.tsv.gz # Link input to output for testing
        '''
}
