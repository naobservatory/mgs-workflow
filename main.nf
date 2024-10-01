include { RUN } from "./workflows/run"
include { RUN2 } from "./workflows/run2"
include { INDEX } from "./workflows/index"
include { BASECALL } from "./workflows/run_ont"

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
    else if (params.mode == "run_ont") {
        RUN_ONT()
    }
}

output {
    directory "${pubDir}"
    mode "copy"
}
