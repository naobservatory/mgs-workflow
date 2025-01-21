// JOIN JSON lengths
process MERGE_JSONS {
    label "biopython"
    input:
        path(input_jsons)
        val(stage_label)
    output:
        path("${stage_label}_length_stats.json")
    shell:
        '''
        merge-json.py -i !{input_jsons} -o !{stage_label}_length_stats.json
        '''
}