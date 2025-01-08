// Summarize Kraken output with Bracken
process BRACKEN {
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
        # Define input/output
        db=!{db_path}
        in=!{report}
        out=!{sample}.bracken
        # Handle gzipped report file
        gzip_in=$(file $(readlink -f ${in}) | grep -c gzip)
        if [[ ${gzip_in} > 0 ]]; then
            echo "Gzipped report - unzipping."
            zcat ${in} > ${in%.gz}
            echo "New input path: ${in%.gz}"
            in=${in%.gz}
        fi
        # Handle empty files and files with no assigned reads
        x=$(wc -l < ${in})
        y=$(grep "unclassified" ${in} | wc -l)
        z=$(cat ${in} | awk '{print $6}' | grep "!{classificationLevel}" | wc -l)
        echo ${x} ${y} ${z}
        if [[ ${x} == "0" && ${y} == "0" ]]; then
            echo "Empty input file - creating empty output."
            touch ${out}
        elif [[ ${x} == "1" && ${y} == "1" ]]; then
            echo "No classified reads in input - creating empty output."
            touch ${out}
        elif [[ ${z} == "0" ]]; then
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
