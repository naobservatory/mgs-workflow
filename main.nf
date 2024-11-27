include { RUN } from "./workflows/run"
include { RUN_VALIDATION } from "./workflows/run_validation"
include { INDEX } from "./workflows/index"

workflow {
    if (params.mode == "index") {
        INDEX()
    } else if (params.mode == "run") {
        RUN()
    } else if (params.mode == "run_validation") {
        RUN_VALIDATION()
    }
}

output {
    "input" {
        path "input"
        tags nextflow_file_class: "publish"
    }
    "logging" {
        path "logging"
        tags nextflow_file_class: "publish"
    }
    "results" {
        path "results"
        tags nextflow_file_class: "publish"
    }
    reads_cleaned {
        path "intermediates/reads/cleaned"
        tags nextflow_file_class: "intermediate"
    }
    reads_raw_viral {
        path "intermediates/reads/raw_viral"
        tags nextflow_file_class: "intermediate"
    }
}
