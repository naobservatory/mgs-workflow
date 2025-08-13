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
include { COUNT_READS_PER_CLADE } from "../modules/local/countReadsPerClade"
include { COPY_FILE_BARE as COPY_VERSION } from "../modules/local/copyFile"
include { COPY_FILE_BARE as COPY_INPUT } from "../modules/local/copyFile"

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
        zero_vv_logs_ch = PREPARE_GROUP_TSVS.out.zero_vv_logs
        // 3. Mark duplicates
        MARK_VIRAL_DUPLICATES(group_ch, params.aln_dup_deviation)
        // Prepare inputs for clade counting and validating taxonomic assignments
        viral_db_path = "${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"
        viral_db = channel.value(viral_db_path)
        dup_ch = MARK_VIRAL_DUPLICATES.out.dup.map{ label, tab, _stats -> [label, tab] }
        // 4. Generate clade counts
        COUNT_READS_PER_CLADE(dup_ch, viral_db)
        // 5. Validate taxonomic assignments
        VALIDATE_VIRAL_ASSIGNMENTS(dup_ch, viral_db,
            params.validation_cluster_identity, 15, params.validation_n_clusters,
            params.ref_dir, params.blast_db_prefix,
            params.blast_perc_id, params.blast_qcov_hsp_perc,
            params.blast_max_rank, params.blast_min_frac,
            params.taxid_artificial)
        // Publish results
        params_str = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
        params_ch = Channel.of(params_str).collectFile(name: "params-downstream.json")
        time_ch = start_time_str.map { it + "\n" }.collectFile(name: "time.txt")
        version_path = file("${projectDir}/pipeline-version.txt")
        version_newpath = version_path.getFileName().toString()
        version_ch = COPY_VERSION(Channel.fromPath(version_path), version_newpath)
        input_newpath = file(params.input_file).getFileName().toString()
        input_file_ch = COPY_INPUT(Channel.fromPath(params.input_file), input_newpath)

    emit:
       input_downstream = params_ch.mix(input_file_ch)
       logging_downstream = time_ch.mix(version_ch)
       intermediates_downstream = VALIDATE_VIRAL_ASSIGNMENTS.out.blast_results
       results_downstream = MARK_VIRAL_DUPLICATES.out.dup.mix(
                                COUNT_READS_PER_CLADE.out.output,
                                VALIDATE_VIRAL_ASSIGNMENTS.out.annotated_hits,
                                zero_vv_logs_ch
                            )
}    
