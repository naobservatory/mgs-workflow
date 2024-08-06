process FASTP {
    label "large"
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
        /* Cleaning not done in CUTADAPT or TRIMMOMATIC:
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
        par="--cut_front --cut_tail --correction --detect_adapter_for_pe --trim_poly_x --cut_mean_quality 25 --average_qual 25 --qualified_quality_phred 20 --verbose --dont_eval_duplication --thread !{task.cpus} --low_complexity_filter"
        # Execute
        fastp ${io} ${par}
        '''
}
