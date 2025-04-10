process EXTRACT_VIRAL_HITS_TO_FASTQ_NOREF {
    label "python"
    label "single"
    input:
        path(tsv) // Viral hits TSV
        val(drop_unpaired) // Boolean
    output:
        path("hits_out.fastq.gz"), emit: output
        path("hits_in.tsv.gz"), emit: input
    shell:
        '''
        extract_viral_hits.py !{drop_unpaired ? "-d" : ""} -i !{tsv} -o hits_out.fastq.gz
        # Link input files for testing
        ln -s !{tsv} hits_in.tsv.gz
        '''
}

process EXTRACT_VIRAL_HITS_TO_FASTQ_NOREF_LABELED {
    label "python"
    label "single"
    input:
        tuple val(sample), path(tsv) // Viral hits TSV
        val(drop_unpaired) // Boolean
    output:
        tuple val(sample), path("${sample}_hits_out.fastq.gz"), emit: output
        tuple val(sample), path("${sample}_hits_in.tsv.gz"), emit: input
    shell:
        '''
        extract_viral_hits.py !{drop_unpaired ? "-d" : ""} -i !{tsv} -o !{sample}_hits_out.fastq.gz
        # Link input files for testing
        ln -s !{tsv} !{sample}_hits_in.tsv.gz
        '''
}
