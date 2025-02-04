// Add a header line to an unheaded TSV file
process HEAD_TSV {
    label "python"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(headers)
        val(label)
    output:
        tuple val(sample), path("${sample}_${label}_head.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_${label}_in.tsv.gz"), emit: input
    shell:
        '''
        head_tsv.py !{tsv} !{headers} !{sample}_!{label}_head.tsv.gz
        # Link input to output for testing
        ln -s !{tsv} !{sample}_!{label}_in.tsv.gz
        '''
}
