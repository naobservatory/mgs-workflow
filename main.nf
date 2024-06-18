include { RUN } from "./workflows/run"

// Configure working and output directories
pubDir  = "${params.base_dir}/output"

workflow {
    RUN()
}

output {
    directory "${pubDir}"
    mode "copy"
}
