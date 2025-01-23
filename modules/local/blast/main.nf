// BLAST paired reads against a downloaded DB and return a single output file
process BLAST_PAIRED_LOCAL {
    label "BLAST"
    label "blast_resources"
    input:
        path(gzipped_reads)
        path(blast_db_dir)
        val(db_prefix)
    output:
        path("hits.blast.gz")
    shell:
        '''
        # Specify parameters
        in1=!{gzipped_reads[0]}
        in2=!{gzipped_reads[1]}
        out=hits.blast
        db=!{blast_db_dir}/!{db_prefix}
        threads=!{task.cpus}
        # Concatenate and decompress inputs
        cat ${in1} ${in2} | zcat > reads.fasta
        # Set up command
        io="-query reads.fasta -out ${out} -db ${db}"
        par="-perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30 -num_threads ${threads}"
        # Run BLAST
        blastn ${io} ${par} -outfmt "6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send"
        # Gzip output
        gzip hits.blast
        '''
}

// BLAST a single input FASTA file against a downloaded DB and return a single output file
process BLAST_LOCAL {
    label "BLAST"
    label "blast_resources"
    input:
        path(gzipped_reads_fasta)
        path(blast_db_dir)
        val(db_prefix)
    output:
        path("hits.blast.gz")
    shell:
        '''
        # Specify BLAST parameters
        in=!{gzipped_reads_fasta}
        out=hits.blast
        db=!{blast_db_dir}/!{db_prefix}
        threads=!{task.cpus}
        # Decompress reads
        zcat ${in} > reads.fasta
        # Set up command
        io="-query reads.fasta -out ${out} -db ${db}"
        par="-perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30 -num_threads ${threads}"
        # Run BLAST
        blastn ${io} ${par} -outfmt "6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send"
        # Gzip output
        gzip ${out}
        '''
}

// BLASTN (streamed version)
process BLASTN_STREAMED {
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
