// Filter genomes to exclude specific patterns in sequence headers
process FILTER_GENOME_FASTA {
    label "seqtk"
    label "single"
    input:
        path(collated_genomes)
        path(patterns_exclude)
    output:
        path("${params.name}.fasta.gz")
    shell:
        '''
        zcat !{collated_genomes} | grep "^>" | grep -vif !{patterns_exclude} | sed 's/>//' > names.txt
        seqtk subseq !{collated_genomes} names.txt | gzip -c > human-viral-genomes-filtered.fasta.gz
        '''
}
