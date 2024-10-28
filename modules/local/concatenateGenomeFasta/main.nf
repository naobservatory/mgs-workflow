// Concatenate downloaded genomes from ncbi-genome-download according to a file of genome IDs
process CONCATENATE_GENOME_FASTA {
    label "single"
    label "seqtk"
    input:
        path(genome_dir)
        path(gid_file)
    output:
        path("genomes.fasta.gz")
    shell:
        '''
        # Get list of matching filenames and save to file
        while read pattern; do find !{genome_dir} -type f -name "*${pattern}*" >> filenames_match.txt; done < !{gid_file}
        # Concatenate
        cat $(cat filenames_match.txt) > genomes.fasta.gz
        '''
}

