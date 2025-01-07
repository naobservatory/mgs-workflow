// 7.3. Perform taxonomic assignment with Kraken2
process KRAKEN {
    label "Kraken2"
    label "kraken_resources"
    input:
        tuple val(sample), path(reads)
        path db_path
    output:
        tuple val(sample), path("${sample}.output"), emit: output
        tuple val(sample), path("${sample}.report"), emit: report
    shell:
        '''
        # Define input/output
        db=!{db_path}
        in=!{reads}
        out=!{sample}.output
        report=!{sample}.report
        io="--output ${out} --report ${report} ${in}"
        # Define parameters
        par="--db ${db} --use-names --report-minimizer-data --threads !{task.cpus}"
        # Run Kraken
        kraken2 ${par} ${io}
        # Make empty output file if no reads
        touch ${out}
        '''
}

// Perform taxonomic assignment with Kraken2 on streamed data
process KRAKEN_STREAMED {
    label "Kraken2"
    label "kraken_resources"
    input:
        tuple val(sample), path(reads)
        path db_path
    output:
        tuple val(sample), path("${sample}.output.gz"), emit: output
        tuple val(sample), path("${sample}.report.gz"), emit: report
        tuple val(sample), path("${sample}_in.fastq.gz"), emit: input
    shell:
        '''
        # Define input/output
        db=!{db_path}
        out=!{sample}.output
        report=!{sample}.report
        # Define parameters
        par="--db ${db} --use-names --report-minimizer-data --threads !{task.cpus} --report ${report}"
        # Run Kraken
        zcat !{reads} | kraken2 ${par} /dev/fd/0 > ${out}
        # Make empty output file if no reads
        touch ${out}
        # Gzip output and report to save space
        gzip ${out}
        gzip ${report}
        # Link input to output for testing
        ln -s !{reads} !{sample}_in.fastq.gz
        '''
}
