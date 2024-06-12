// Subsample reads with seqtk
process SUBSET_READS_PAIRED {
    label "seqtk"
    label "single"
    input:
        tuple val(sample), path(reads)
        val readFraction
    output:
        tuple val(sample), path("${sample}_subset_{1,2}.fastq.gz")
    shell:
        '''
        # Define input/output
        in1=!{reads[0]}
        in2=!{reads[1]}
        out1=!{sample}_subset_1.fastq.gz
        out2=!{sample}_subset_2.fastq.gz
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
