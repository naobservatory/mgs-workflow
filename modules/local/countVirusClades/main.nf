// Extract table of clade counts from virus reads
process COUNT_VIRUS_CLADES {
    label "R"
    label "large"
    input:
        path(virus_hits_tsv)
        path(virus_db)
    output:
        path("virus_clade_counts.tsv.gz")
    shell:
        '''
        count-viral-taxa.R --reads !{virus_hits_tsv} --taxa !{virus_db} --output virus_clade_counts.tsv.gz
        '''
}
