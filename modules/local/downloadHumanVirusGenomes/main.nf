// Get viral genome sequences based on expanded list of taxids
process DOWNLOAD_HUMAN_VIRUS_GENOMES {
    label "biopython"
    label "max"
    input:
        path(hv_taxids)
    output:
        path("ncbi_fetch_metadata.txt"), emit: metadata
        path("genbank_genomes"), emit: genomes
    shell:
        '''
        download-ncbi-genomes.py !{hv_taxids} --threads !{task.cpus} --output genbank_genomes --metadata ncbi_fetch_metadata.txt
        '''
}
