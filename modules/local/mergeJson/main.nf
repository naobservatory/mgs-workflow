// JOIN JSON lengths
process GET_READ_LENGTHS {
    label "biopython"
    input:
        path(input_jsons)
    output:
        path("read_lengths.json")
    shell:
        '''
        merge-json.py -i !{input_jsons} -o "read_lengths.json"
        '''
}