process ATRIA {
    label "atria"
    label "large"
    input:
        // reads is a list of two files: forward/reverse reads
        tuple val(sample), path(reads)
        val(adapters_ch)
    output:
        tuple val(sample), path("*{1,2}.atria.fastq.gz"), emit: reads
    shell:
        '''
        par="--adapter1 !{adapters_ch.join(' ')} --adapter2 !{adapters_ch.join(' ')} --kmer-tolerance 2 --kmer-n-match 9 --trim-score-pe 10.0 --quality-kmer 4 --quality-score 15 --length-range 20:999999 --no-consensus"
        in="--read1 !{reads[0]} --read2 !{reads[1]}"
        atria ${par} ${in}
        '''
}
