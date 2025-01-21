// Quantify read lengths
process GET_READ_LENGTHS {
    label "biopython"
    input:
        tuple val(sample), path(input_fastq)
        val(stage)
    output:
        path("${stage}_${sample}_read_lengths.json")
    shell:
        '''
        get-read-lengths.py -i !{input_fastq} -stage !{stage} -sample !{sample} -o "!{stage}_!{sample}_read_lengths.json"
        '''
}