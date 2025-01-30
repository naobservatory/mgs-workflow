// Subsample single or interleaved FASTQ based on read IDs
// NB: Currently assumes gzipped input
process SUBSET_FASTQ {
    label "seqtk"
    label "single"
    input:
        tuple val(sample), path(fastq)
        val readFraction
        val randomSeed
    output:
        tuple val(sample), path("${sample}_reads_subset.fastq.gz"), emit: output
        tuple val(sample), path("${sample}_reads_in.fastq.gz"), emit: input
    shell:
        '''
        rf=!{readFraction}
        out=!{sample}_reads_subset.fastq.gz
        in=!{sample}_reads_in.fastq.gz
        # Check the read fraction is between 0 and 1 inclusive
        if ! awk -v f="${rf}" 'BEGIN {if (f < 0 || f > 1) exit 1}'; then
            echo "ERROR: Read fraction must be between 0 and 1 (got ${rf})."
            exit 1
        fi
        # Handle special cases where fraction equals 0 or 1
        if awk -v f="${rf}" 'BEGIN {if (f == 1) exit 0; else exit 1}'; then
            echo "Read fraction equals 1. Copying input to output without sampling."
            cp !{fastq} ${out}
            ln -s !{fastq} ${in}
            exit 0
        elif awk -v f="${rf}" 'BEGIN {if (f == 0) exit 0; else exit 1}'; then
            echo "Read fraction equals 0. Creating empty output file."
            touch ${out%.gz}
            gzip ${out%.gz}
            ln -s !{fastq} ${in}
            exit 0
        fi
        # Get unique read IDs from input FASTQ file
        zcat !{fastq} | awk 'NR % 4 == 1 {
            sub(/^@/, "")
            split($0, fields, " ")
            print fields[1]
        }' | sort -u > all_ids.txt
        echo "Input IDs $(cat all_ids.txt | wc -l)"
        # Subsample IDs
        frac=$(echo | awk -v f=${rf} '{ ff = f + 0; print ff }')
        echo "Target read fraction: ${frac}"
        rseed=!{randomSeed == "" ? "\\$RANDOM" : randomSeed}
        echo "Random seed: ${rseed}"
        awk -v f=${rf} -v s=${rseed} '
            BEGIN {frac=f+0;srand(s);} {if (rand() <= frac) {print}}
        ' all_ids.txt > sample_ids.txt
        echo "Output IDs $(cat sample_ids.txt | wc -l)"
        # Subsample sequences based on IDs
        seqtk subseq !{fastq} sample_ids.txt | gzip -c > ${out}
        # Link input to output for testing
        ln -s !{fastq} ${in}
        '''
}
