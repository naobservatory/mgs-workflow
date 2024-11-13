// Download entire viral Genbank DB
process DOWNLOAD_VIRAL_GENBANK {
    label "biopython"
    label "max"
    output:
        path("genbank_metadata.txt"), emit: metadata
        path("genbank_genomes"), emit: genomes
    shell:
        '''
        par="--section genbank --formats fasta --flat-output --verbose --parallel !{task.cpus}"
        io="--output-folder genbank_genomes --metadata-table genbank_metadata.txt"
        ncbi-genome-download ${par} ${io} viral
        '''
}
