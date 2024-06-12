// Summarize Kraken output with Bracken
process BRACKEN {
    label "Bracken"
    label "single"
    input:
        tuple val(sample), path(report)
        path db_path
        val classificationLevel
    output:
        tuple val(sample), path("${sample}.bracken")
    shell:
        '''
        # Define input/output
        db=!{db_path}
        in=!{report}
        out=!{sample}.bracken
        io="-d ${db} -i ${in} -o ${out}"
        # Define parameters
        par="-l !{classificationLevel}"
        # Handle files with no assigned reads
        x=$(wc -l < ${in})
        y=$(grep "unclassified" ${in} | wc -l)
        echo ${x} ${y}
        if [[ ${x} == "1" && ${y} == "1" ]]; then
            touch ${out}
        else
            # Run Bracken
            bracken ${io} ${par}
        fi
        '''
}
