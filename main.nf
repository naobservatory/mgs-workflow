include { RUN } from "./workflows/run"
include { RUN2 } from "./workflows/run2"
include { INDEX } from "./workflows/index"
include { BASECALL } from "./workflows/basecall"

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
    else if (params.mode == "basecall") {
        BASECALL()
    }
}

output {
    directory "${pubDir}"
    mode "copy"
}
