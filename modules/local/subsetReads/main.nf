// Subsample reads with seqtk
process SUBSET_READS_PAIRED {
    label "seqtk"
    label "single"
    input:
        tuple val(sample), path(reads)
        val readFraction
    output:
        tuple val(sample), path("${sample}_subset_{1,2}.${params.suffix}.gz")
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        out1=!{sample}_subset_1.!{params.suffix}.gz
        out2=!{sample}_subset_2.!{params.suffix}.gz
        # Count reads for validation
        echo "Input reads: $(zcat ${in1} | wc -l | awk '{ print $1/4 }')"
        # Carry out subsetting
        seed=${RANDOM}
        seqtk sample -s ${seed} ${in1} !{readFraction} | gzip -c > ${out1}
        seqtk sample -s ${seed} ${in2} !{readFraction} | gzip -c > ${out2}
        # Count reads for validation
        echo "Output reads: $(zcat ${out1} | wc -l | awk '{ print $1/4 }')"
        '''
}

// Subsample reads with seqtk (no sample name)
process SUBSET_READS_PAIRED_MERGED {
    label "seqtk"
    label "single"
    input:
        path(reads)
        val readFraction
    output:
        path("reads_subset_{1,2}.${params.suffix}.gz")
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        out1=reads_subset_1.!{params.suffix}.gz
        out2=reads_subset_2.!{params.suffix}.gz
        # Count reads for validation
        echo "Input reads: $(zcat ${in1} | wc -l | awk '{ print $1/4 }')"
        # Carry out subsetting
        seed=${RANDOM}
        seqtk sample -s ${seed} ${in1} !{readFraction} | gzip -c > ${out1}
        seqtk sample -s ${seed} ${in2} !{readFraction} | gzip -c > ${out2}
        # Count reads for validation
        echo "Output reads: $(zcat ${out1} | wc -l | awk '{ print $1/4 }')"
        '''
}
