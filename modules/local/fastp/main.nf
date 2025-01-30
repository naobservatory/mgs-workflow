// Run FASTP on streamed data (either single-end or interleaved)
process FASTP {
    label "small"
    label "fastp"
    input:
        tuple val(sample), path(reads)
        path(adapters)
        val(interleaved)
    output:
        tuple val(sample), path("${sample}_fastp.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_fastp_failed.fastq.gz"), emit: failed
        tuple val(sample), path("${sample}_fastp.{json,html}"), emit: log
        tuple val(sample), path("${sample}_fastp_in.fastq.gz"), emit: input
    shell:
        /* Cleaning not done in CUTADAPT:
        * Higher quality threshold for sliding window trimming;
        * Removing poly-X tails;
        * Automatic adapter detection;
        * Base correction in overlapping paired-end reads;
        * Filter low complexity reads.
        */
        '''
        # Define paths and parameters
        op=!{sample}_fastp.fastq.gz
        of=!{sample}_fastp_failed.fastq.gz
        oj=!{sample}_fastp.json
        oh=!{sample}_fastp.html
        ad=!{adapters}
        io="--failed_out ${of} --html ${oh} --json ${oj} --adapter_fasta ${ad} --stdin --stdout !{interleaved ? '--interleaved_in' : ''}"
        par="--cut_front --cut_tail --correction --detect_adapter_for_pe --trim_poly_x --cut_mean_quality 20 --average_qual 20 --qualified_quality_phred 20 --verbose --dont_eval_duplication --thread !{task.cpus} --low_complexity_filter"
        # Execute
        zcat !{reads} | fastp ${io} ${par} | gzip -c > ${op}
        # Link input to output for testing
        ln -s !{reads} !{sample}_fastp_in.fastq.gz
        '''
}

// Run FASTP for adapter trimming but don't trim for quality
process FASTP_NOTRIM {
    label "max"
    label "fastp"
    input:
        // reads is a list of two files: forward/reverse reads
        tuple val(sample), path(reads)
        path(adapters)
    output:
        tuple val(sample), path("${sample}_fastp_{1,2}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_fastp_failed.fastq.gz"), emit: failed
        tuple val(sample), path("${sample}_fastp.{json,html}"), emit: log
    shell:
        /* Cleaning not done in CUTADAPT:
        * Higher quality threshold for sliding window trimming;
        * Removing poly-X tails;
        * Automatic adapter detection;
        * Base correction in overlapping paired-end reads;
        * Filter low complexity reads.
        */
        '''
        # Define paths and subcommands
        o1=!{sample}_fastp_1.fastq.gz
        o2=!{sample}_fastp_2.fastq.gz
        of=!{sample}_fastp_failed.fastq.gz
        oj=!{sample}_fastp.json
        oh=!{sample}_fastp.html
        ad=!{adapters}
        io="--in1 !{reads[0]} --in2 !{reads[1]} --out1 ${o1} --out2 ${o2} --failed_out ${of} --html ${oh} --json ${oj} --adapter_fasta ${ad}"
        par="--detect_adapter_for_pe --disable_quality_filtering --disable_length_filtering --verbose --dont_eval_duplication --thread !{task.cpus}"
        # Execute
        fastp ${io} ${par}
        '''
}

