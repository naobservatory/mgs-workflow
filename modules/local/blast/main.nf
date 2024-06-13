// BLAST paired reads against NT and return a single output file
process BLAST_PAIRED_NT {
    label "BLAST"
    cpus "${params.cpus}"
    memory "${params.mem}"
    input:
        path(reads)
        path(blast_nt_dir)
    output:
        path("hits.blast.gz")
    shell:
        '''
        # Specify parameters
        in1=!{reads[0]}
        in2=!{reads[1]}
        out=hits.blast
        db=!{blast_nt_dir}/nt
        threads=!{task.cpus}
        # Concatenate inputs
        cat ${in1} ${in2} > reads.fasta
        # Set up command
        io="-query reads.fasta -out ${out} -db ${db}"
        par="-perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30 -num_threads ${threads}"
        # Run BLAST
        blastn ${io} ${par} -outfmt "6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send"
        # Gzip output
        gzip hits.blast
        '''
}
