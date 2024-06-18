include { RUN } from "./workflows/run"
include { INDEX } from "./workflows/index"

// Configure working and output directories
pubDir  = "${params.base_dir}/output"

workflow {
    if (params.mode == "index") {
        INDEX()
    } else if (params.mode == "run") {
        RUN()
    }
}

output {
    directory "${pubDir}"
    mode "copy"
}
