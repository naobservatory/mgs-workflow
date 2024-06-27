// Label Bracken files with sample IDs
process LABEL_BRACKEN_REPORTS {
    label "tidyverse"
    label "single"
    input:
        tuple val(sample), path(bracken_output)
    output:
        path("${sample}_labeled.bracken")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        tab <- read_tsv("!{bracken_output}", show_col_types = FALSE) %>%
            mutate(sample="!{sample}")
        write_tsv(tab, "!{sample}_labeled.bracken")
        '''
}
