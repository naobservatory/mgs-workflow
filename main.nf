include { RUN } from "./workflows/run"

workflow {
    RUN()
}

output {
    directory "${params.pub_dir}"
    mode "copy"
}
