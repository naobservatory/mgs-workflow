include { RUN } from "./workflows/run"
include { RUN2 } from "./workflows/run2"
include { RUN_VALIDATION } from "./workflows/run_validation"
include { INDEX } from "./workflows/index"

// Configure working and output directories
pubDir  = "${params.base_dir}/output"

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
    directory "${pubDir}"
    mode "copy"
}
