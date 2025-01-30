// BLASTN (streamed version)
process BLASTN {
    label "BLAST"
    label "blast_resources"
    input:
        tuple val(sample), path(fasta) // Gzipped interleaved or single-end FASTA
        path(blast_db_dir)
        val(db_prefix)
    output:
        tuple val(sample), path("${sample}_hits.blast.gz"), emit: output
        tuple val(sample), path("${sample}_in.fasta.gz"), emit: input
    shell:
        '''
        # Set up command
        io="-db !{blast_db_dir}/!{db_prefix}"
        par="-perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30 -num_threads !{task.cpus}"
        fmt="6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send"
        # Run BLAST
        zcat !{fasta} | blastn ${io} ${par} -outfmt "${fmt}" \\
            | gzip > !{sample}_hits.blast.gz
        # Link input to output for testing
        ln -s !{fasta} !{sample}_in.fasta.gz
        '''
}
