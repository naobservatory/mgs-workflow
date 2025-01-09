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

// Summarize Kraken output with Bracken (updated version)
process BRACKEN2 {
    label "bracken_plus_utils"
    label "single"
    input:
        tuple val(sample), path(report)
        path db_path
        val classificationLevel
    output:
        tuple val(sample), path("${sample}.bracken.gz")
    shell:
        '''
        #set -euox pipefail
        # Define input/output
        db=!{db_path}
        in=!{report}
        out=!{sample}.bracken
        # Handle gzipped report file
        in_status=$(file -b $(readlink -f ${in}))
        echo "Input file status: ${in_status}"
        gzip_in=$(echo "${in_status}" | grep -c "gzip" || true)
        echo "Gzip count in file status: ${gzip_in}"
        if [[ ${gzip_in} -gt 0 ]]; then
            echo "Gzipped report - unzipping."
            zcat ${in} > ${in%.gz}
            echo "New input path: ${in%.gz}"
            in=${in%.gz}
        else
            echo "Plain input - no unzipping needed."
        fi
        # Handle empty files and files with no assigned reads
        x=$(wc -l < ${in})
        y=$(grep -c "unclassified" ${in} || true)
        z=$(cat "${in}" | awk '{print $6}' | grep -c "!{classificationLevel}" || true)
        echo "Number of lines (x): ${x}"
        echo "Number of 'unclassified' lines (y): ${y}"
        echo "Number of lines with classification level !{classificationLevel} (z): ${z}"
        if [[ ${x} -eq "0" && ${y} -eq "0" ]]; then
            echo "Empty input file - creating empty output."
            touch ${out}
        elif [[ ${x} -eq "1" && ${y} -eq "1" ]]; then
            echo "No classified reads in input - creating empty output."
            touch ${out}
        elif [[ ${z} -eq "0" ]]; then
            echo "No reads classified at desired level in input - creating empty output."
            touch ${out}
        else
            # Run Bracken
            io="-d ${db} -i ${in} -o ${out}"
            par="-l !{classificationLevel}"
            echo "Input okay - running Bracken."
            bracken ${io} ${par}
        fi
        # Gzip output
        gzip -c ${out} > ${out}.gz
        '''
}
