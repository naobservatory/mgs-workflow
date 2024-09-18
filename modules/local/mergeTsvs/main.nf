// Combine multiple TSVs with identical headers into a single output file
process MERGE_TSVS {
    label "tidyverse"
    cpus 1
    memory "16.GB"
    input:
        path(tsvs)
    output:
        path("${params.name}.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        in_paths <- str_split("!{tsvs}", " ")[[1]]
        print(in_paths)
        tabs <- lapply(in_paths, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        for (t in tabs) print(dim(t))
        tab_out <- bind_rows(tabs)
        print(dim(tab_out))
        sapply(tabs, nrow) %>% sum %>% print
        write_tsv(tab_out, "!{params.name}.tsv.gz")
        '''
}
