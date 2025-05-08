/***********************************************************
| WORKFLOW: DOWNSTREAM ANALYSIS OF PRIMARY WORKFLOW OUTPUT |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { LOAD_DOWNSTREAM_DATA } from "../subworkflows/local/loadDownstreamData"
include { PREPARE_GROUP_TSVS } from "../subworkflows/local/prepareGroupTsvs"
include { MARK_VIRAL_DUPLICATES } from "../subworkflows/local/markViralDuplicates"
include { VALIDATE_VIRAL_ASSIGNMENTS } from "../subworkflows/local/validateViralAssignments"

nextflow.preview.output = true

/*****************
| MAIN WORKFLOWS |
*****************/

workflow DOWNSTREAM {
    main:
        // 1. Prepare channel from input TSV file
        LOAD_DOWNSTREAM_DATA(params.input_file)
        input_ch = LOAD_DOWNSTREAM_DATA.out.input
        start_time_str = LOAD_DOWNSTREAM_DATA.out.start_time_str
        // 2. Add group information, partition into per-group TSVs
        PREPARE_GROUP_TSVS(input_ch)
        group_ch = PREPARE_GROUP_TSVS.out.groups
        // 3. Mark duplicates
        MARK_VIRAL_DUPLICATES(group_ch, params.aln_dup_deviation)
        // 4. Validate taxonomic assignments
        viral_db_path = "${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"
        viral_db = Channel.of(viral_db_path)
        dup_ch = MARK_VIRAL_DUPLICATES.out.dup.map{ label, tab, stats -> [label, tab] }
        //VALIDATE_VIRAL_ASSIGNMENTS(dup_ch, viral_db,
        //    params.validation_cluster_identity, 15, params.validation_n_clusters, params.single_end)
        // Publish results
        params_str = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
        params_ch = Channel.of(params_str).collectFile(name: "params-downstream.json")
        time_ch = start_time_str.map { it + "\n" }.collectFile(name: "time.txt")
        version_ch = Channel.fromPath("${projectDir}/pipeline-version.txt")
    publish:
        // Saved inputs
        Channel.fromPath(params.input_file) >> "input_downstream"
        params_ch >> "input_downstream"
        time_ch >> "logging_downstream"
        version_ch >> "logging_downstream"
        // Duplicate results
        MARK_VIRAL_DUPLICATES.out.dup >> "results_downstream"
}
