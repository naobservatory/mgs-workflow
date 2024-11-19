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
    } else if (params.mode == "test") {
        INDEX()
        RUN()
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
