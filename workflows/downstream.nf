/***********************************************************
| WORKFLOW: DOWNSTREAM ANALYSIS OF PRIMARY WORKFLOW OUTPUT |
***********************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { LOAD_DOWNSTREAM_DATA } from "../subworkflows/local/loadDownstreamData"
include { PREPARE_GROUP_TSVS } from "../subworkflows/local/prepareGroupTsvs"

nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

workflow DOWNSTREAM {

    // 1. Prepare channel from input TSV file
    LOAD_DOWNSTREAM_DATA(params.input_file)
    input_ch = LOAD_DOWNSTREAM_DATA.out.input
    start_time_str = LOAD_DOWNSTREAM_DATA.out.start_time_str

    // 2. Add group information, partition into per-group TSVs
    PREPARE_GROUP_TSVS(input_ch)
    group_ch = PREPARE_GROUP_TSVS.out.groups

//    // 3. Mark duplicates
//    MARK_DUPLICATES(group_ch)

    // Publish results
    params_str = JsonOutput.prettyPrint(JsonOutput.toJson(params))
    params_ch = Channel.of(params_str).collectFile(name: "downstream-params.json")
    time_ch = start_time_str.map { it + "\n" }.collectFile(name: "time.txt")
    version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    publish:
        // Saved inputs
        Channel.fromPath(params.input_file) >> "input"
        params_ch >> "input"
        time_ch >> "logging"
        version_ch >> "logging"
        // Duplicate results
//        MARK_DUPLICATES.out.tsv >> "results_downstream"
//        MARK_DUPLICATES.out.mapping >> "results_downstream"
}
