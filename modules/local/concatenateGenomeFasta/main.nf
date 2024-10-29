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
        # Diagnostics
        echo "Genome directory contains" $(ls !{genome_dir} | wc -l) "files, beginning with:"
        ls -1 !{genome_dir} | head
        echo "GID file contains" $(cat !{gid_file} | wc -l) "patterns, beginning with:"
        head !{gid_file}
        # Get list of matching filenames and save to file
        while read pattern; do 
            ls -1 !{genome_dir}/*${pattern}* >> filenames_match.txt
            #find !{genome_dir} -type f -name "*${pattern}*" >> filenames_match.txt
        done < !{gid_file}
        echo "Temporary file contains" $(cat filenames_match.txt | wc -l) "filepaths, beginning with:"
        head filenames_match.txt
        # Concatenate
        if [[ -s filenames_match.txt ]]; then
            cat $(cat filenames_match.txt) > genomes.fasta.gz
        else
            echo "No matching files found!"
            exit 1
        fi
        echo "Output file contains" $(zcat genomes.fasta.gz | grep "^>" | wc -l) "sequences."
        '''
}
