// Perform taxonomic assignment with Kraken2 on streamed data
process KRAKEN {
    label "Kraken2"
    label "kraken_resources"
    input:
        tuple val(sample), path(reads)
        val db_path
    output:
        tuple val(sample), path("${sample}.output.gz"), emit: output
        tuple val(sample), path("${sample}.report.gz"), emit: report
        tuple val(sample), path("${sample}_in.fastq.gz"), emit: input
    shell:
        '''
        download-db.sh !{db_path}
        #cp -r /fsx/20250404_index_output/results/kraken_db /scratch/kraken_db && touch /scratch/kraken_db/download_complete.txt
        # Define input/output
        out=!{sample}.output
        report=!{sample}.report
        # Define parameters
        db_name=\$(basename "!{db_path}")
        par="--db /scratch/\${db_name} --use-names --report-minimizer-data --threads !{task.cpus} --report ${report} --memory-mapping"
        # Run Kraken
        zcat !{reads} | kraken2 ${par} /dev/fd/0 > ${out}
        # Make empty output files if needed
        touch ${out}
        touch ${report}
        # Gzip output and report to save space
        gzip ${out}
        gzip ${report}
        # Link input to output for testing
        ln -s !{reads} !{sample}_in.fastq.gz
        '''
}
