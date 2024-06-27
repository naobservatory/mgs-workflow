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
        # Handle empty files and files with no assigned reads
        x=$(wc -l < ${in})
        y=$(grep "unclassified" ${in} | wc -l)
        z=$(cat ${in} | awk '{print $6}' | grep "!{classificationLevel}" | wc -l)
        echo ${x} ${y}
        if [[ ${x} == "0" && ${y} == "0" ]]; then
            touch ${out}
        elif [[ ${x} == "1" && ${y} == "1" ]]; then
            touch ${out}
        elif [[ ${z} == "0" ]]; then
            touch ${out}
        else
            # Run Bracken
            bracken ${io} ${par}
        fi
        '''
}
