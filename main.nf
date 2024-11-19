include { RUN } from "./workflows/run"
include { RUN_VALIDATION } from "./workflows/run_validation"
include { INDEX } from "./workflows/index"
include { RUN_DEV_SE } from "./workflows/run_dev_se"

workflow {
    if (params.mode == "index") {
        INDEX()
    } else if (params.mode == "run") {
        RUN()
    } else if (params.mode == "run_validation") {
        RUN_VALIDATION()
    } else if (params.mode == "run_dev_se") {
        RUN_DEV_SE()
    }
}

output {
    "input" {
        path "input"
    }
    "logging" {
        path "logging"
    }
    "results" {
        path "results"
    }
}
