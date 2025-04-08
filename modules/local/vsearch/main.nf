// Cluster merged FASTQ sequences with VSEARCH
process VSEARCH_CLUSTER {
    label "large"
    label "vsearch"
    input:
        tuple val(sample), path(reads) // Single-end or merged reads
        val(identity_threshold) // Minimum required identity (0.0-1.0) required for two sequences to cluster together
        val(identity_method) // Method for calculating identity (see VSEARCH documentation)
        val(min_seq_length) // Minimum sequencing length required by VSEARCH
    output:
        tuple val(sample), path("${sample}_vsearch_reps.fasta.gz"), emit: reps
        tuple val(sample), path("${sample}_vsearch_summary.tsv.gz"), emit: summary
        tuple val(sample), path("${sample}_vsearch_log.txt"), emit: log
        tuple val(sample), path("input_${reads}"), emit: input
    shell:
        '''
        # Define paths and parameters
        or=!{sample}_vsearch_reps.fasta
        os=!{sample}_vsearch_summary.tsv
        log=!{sample}_vsearch_log.txt
        io="--log ${log} --centroids ${or} --uc ${os} --cluster_fast !{reads}"
        par="--threads !{task.cpus} --id !{identity_threshold} --iddef !{identity_method} --minseqlength !{min_seq_length}"
        # Add decompression if necessary
        par="${par}!{reads.endsWith(".gz") ? ' --gzip_decompress' : ''}"
        # Execute
        vsearch ${par} ${io}
        # Gzip outputs
        gzip -c ${or} > ${or}.gz
        gzip -c ${os} > ${os}.gz
        # Link input to output for testing
        ln -s !{reads} input_!{reads}
        '''
}
