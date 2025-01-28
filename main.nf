include { RUN } from "./workflows/run"
include { RUN_VALIDATION } from "./workflows/run_validation"
include { INDEX } from "./workflows/index"
include { RUN_DEV_SE } from "./workflows/run_dev_se"
include { RUN_STREAMED } from "./workflows/run_streamed"
include { RUN_VALIDATION_STREAMED } from "./workflows/run_validation_streamed"

workflow {
    if (params.mode == "index") {
        INDEX()
    } else if (params.mode == "run") {
        RUN()
    } else if (params.mode == "run_validation") {
        RUN_VALIDATION()
    } else if (params.mode == "run_dev_se") {
        RUN_DEV_SE()
    } else if (params.mode == "run_streamed") {
        RUN_STREAMED()
    } else if (params.mode == "run_validation_streamed") {
        RUN_VALIDATION_STREAMED()
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
    intermediates {
        path "intermediates"
        tags nextflow_file_class: "intermediate", "nextflow.io/temporary": "false"
    }
}
