include { RUN } from "./workflows/run"
include { RUN_VALIDATION } from "./workflows/run_validation"
include { INDEX } from "./workflows/index"
include { RUN_DEV_SE } from "./workflows/run_dev_se"
include { DOWNSTREAM } from "./workflows/downstream"

workflow {
    if (params.mode == "index") {
        INDEX()
    } else if (params.mode == "run") {
        RUN()
    } else if (params.mode == "run_validation") {
        RUN_VALIDATION()
    } else if (params.mode == "run_dev_se") {
        RUN_DEV_SE()
    } else if (params.mode == "downstream") {
        DOWNSTREAM()
    }
}

output {
    "input" {
        path "input"
        tags nextflow_file_class: "publish", "nextflow.io/temporary": "false"
    }
    "logging" {
        path "logging"
        tags nextflow_file_class: "publish", "nextflow.io/temporary": "false"
    }
    "results" {
        path "results"
        tags nextflow_file_class: "publish", "nextflow.io/temporary": "false"
    }
    reads_raw_viral {
        path "intermediates/reads/raw_viral"
        tags nextflow_file_class: "intermediate", "nextflow.io/temporary": "false"
    }
    reads_trimmed_viral {
        path "intermediates/reads/trimmed_viral"
        tags nextflow_file_class: "intermediate", "nextflow.io/temporary": "false"
    }
    intermediates {
        path "intermediates"
        tags nextflow_file_class: "intermediate", "nextflow.io/temporary": "false"
    }
    "results_downstream" {
        path "results_downstream"
        tags nextflow_file_class: "publish", "nextflow.io/temporary": "false"
    }
}
