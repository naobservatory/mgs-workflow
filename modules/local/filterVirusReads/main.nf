// Filter virus reads by alignment score and assignment status (streamed version)
process FILTER_VIRUS_READS {
    label "python"
    label "single"
    input:
        path(hits_tsv)
        val(score_threshold) // Length-normalized Bowtie2 alignment score
    output:
        path("virus_hits_filtered.tsv.gz"), emit: output
        path("virus_hits_in.tsv.gz"), emit: input
    shell:
        '''
        filter_virus_reads.py !{hits_tsv} !{score_threshold} virus_hits_filtered.tsv.gz
        ln -s !{hits_tsv} virus_hits_in.tsv.gz
        '''
}
