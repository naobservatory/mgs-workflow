// Concatenate split files from same sample together
process CONCAT_GZIPPED {
    label "single"
    label "seqtk"
    input:
        path(raw_files_directory)
        tuple val(sample), val(libraries)
    output:
        tuple val(sample), path("${sample}_{1,2}.fastq.gz"), emit: reads
    shell:
        '''
        # Preamble
        read_dir=!{raw_files_directory}
        echo Raw files directory: $read_dir
        echo Sample: !{sample}
        echo Libraries: !{libraries.join(" ")}
        # Get file paths from library IDs
        r1=""
        r2=""
        for l in !{libraries.join(" ")}; do
            L1=$(ls ${read_dir}/*${l}*!{params.r1_suffix}.fastq.gz)
            L2=$(ls ${read_dir}/*${l}*!{params.r2_suffix}.fastq.gz)
            r1="${r1} ${L1}"
            r2="${r2} ${L2}"
            done
        ln1=$(wc -w <<< ${r1})
        ln2=$(wc -w <<< ${r2})
        echo Forward read files: ${ln1}
        echo Reverse read files: ${ln2}
        if [[ ${ln1} != ${ln2} ]]; then
            >&2 echo "Error: uneven numbers of files for forward and reverse reads."
        fi
        if [[ ${ln1} == 0 ]]; then
            >&2 echo "Error: No read files specified"
        fi
        # Concatenate or copy
        out1=!{sample}_1.fastq.gz
        out2=!{sample}_2.fastq.gz
        echo Read 1 files to concatenate: ${r1}
        echo Read 2 files to concatenate: ${r2}
        if [[ ${ln1} == 1 ]]; then
            # Make copies
            echo "Only one file per read pair; copying."
            cp ${r1} ${out1}
            cp ${r2} ${out2}
            # Test copies and fail if not identical
            cmp ${r1} ${out1} > diff1.txt
            cmp ${r2} ${out2} > diff2.txt
            if [[ -s diff1.txt ]]; then
                >&2 echo "Error: $(cat diff1)"
                echo "Input file lengths: $(zcat ${r1} | wc -l) $(zcat ${r2} | wc -l)
                echo "Output file lengths: $(zcat ${out1} | wc -l) $(zcat ${out2} | wc -l)
                exit 1
            elif [[ -s diff2.txt ]]; then
                >&2 echo "Error: $(cat diff2)"
                echo "Input file lengths: $(zcat ${r1} | wc -l) $(zcat ${r2} | wc -l)
                echo "Output file lengths: $(zcat ${out1} | wc -l) $(zcat ${out2} | wc -l)
                exit 1
            fi
        else
            # Concatenate
            # TODO: Add error checking and handling
            cat ${r1} > ${out1}
            cat ${r2} > ${out2}
        fi
        '''
}
