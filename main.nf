include { RUN } from "./workflows/run"
include { RUN2 } from "./workflows/run2"
include { INDEX } from "./workflows/index"

// Configure working and output directories
pubDir  = "${params.base_dir}/output"

workflow {
    if (params.mode == "index") {
        INDEX()
    } else if (params.mode == "run") {
        RUN()
    } else if (params.mode == "run2") {
        RUN2()
    }
}

output {
    directory "${pubDir}"
    mode "copy"
}
