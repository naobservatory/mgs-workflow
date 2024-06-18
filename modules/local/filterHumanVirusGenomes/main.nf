// Filter HV genomes to exclude specific patterns (e.g. transgenic, mutated, unverified, or contaminated sequences)
process FILTER_HUMAN_VIRUS_GENOMES {
    label "seqtk"
    label "single"
    input:
        path(collated_genomes)
        path(patterns_exclude)
    output:
        path("human-viral-genomes-filtered.fasta.gz")
    shell:
        '''
        zcat !{collated_genomes} | grep ^> | grep -vif !{patterns_exclude} | sed 's/>//' > names.txt
        seqtk subseq !{collated_genomes} names.txt | gzip -c > human-viral-genomes-filtered.fasta.gz
        '''
}
