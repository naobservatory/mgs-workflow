// 7.3. Perform taxonomic assignment with Kraken2
process KRAKEN {
    label "Kraken2"
    cpus 16
    memory "${mem}"
    input:
        tuple val(sample), path(reads)
        path db_path
        val mem
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
