// Filter virus reads by alignment score and assignment status (streamed version)
process FILTER_VIRUS_READS {
    label "python"
    label "single"
    input:
        path(hits_tsv)
        val(label) // Output filename (needed for publishing)
    output:
        path("${label}.tsv.gz"), emit: output
        path("input_${hits_tsv}"), emit: input
    shell:
        '''
        filter_virus_reads.py !{hits_tsv} !{label}.tsv.gz
        ln -s !{hits_tsv} input_!{hits_tsv}
        '''
}
