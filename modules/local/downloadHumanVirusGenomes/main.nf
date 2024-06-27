// Get viral genome sequences based on expanded list of taxids
process DOWNLOAD_HUMAN_VIRUS_GENOMES {
    label "biopython"
    label "single"
    input:
        path(hv_taxids)
    output:
        path("ncbi_fetch_metadata.txt"), emit: metadata
        path("genbank_genomes"), emit: genomes
    shell:
        '''
        io="--taxids !{hv_taxids} --metadata-table ncbi_fetch_metadata.txt -o genbank_genomes"
        par="--section genbank --formats fasta --flat-output"
        ncbi-genome-download ${par} ${io} viral
        '''
}
