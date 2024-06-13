// Extract table of clade counts from HV reads
process COUNT_HV_CLADES {
    label "R"
    label "large"
    publishDir "${pubDir}", mode: "copy", overwrite: "true"
    input:
        path(hv_hits_tsv)
        path(viral_taxa)
    output:
        path("hv_clade_counts.tsv.gz")
    shell:
        '''
        count-viral-taxa.R --reads !{hv_hits_tsv} --taxa !{viral_taxa} --output hv_clade_counts.tsv.gz
        '''
}
