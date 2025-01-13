process COUNT_READS {
    label "countReads"
    input:
        tuple val(sample), path(reads)
    output:
        path("${sample}_read_count.txt")
    shell:
        '''
        # Count reads in file
        COUNT=$(bioawk -c fastx 'END{print NR}' !{reads[0]})
        echo "!{sample},$COUNT" > !{sample}_read_count.txt
        '''
}