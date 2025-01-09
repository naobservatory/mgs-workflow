// Combine multiple TSVs with identical headers into a single output file
process CONCATENATE_TSVS {
    label "tidyverse"
    label "single_cpu_16GB_memory"
    input:
        path(tsvs)
        val(name)
    output:
        path("${name}.tsv.gz")
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
        write_tsv(tab_out, "!{name}.tsv.gz")
        '''
}

// Concatenate multiple TSVs (streamed version with Python)
process CONCATENATE_TSVS_STREAMED {
    label "biopython"
    label "single"
    input:
        path(tsvs)
        val(name)
    output:
        path("${name}.tsv.gz"), emit: output
        path("${name}_in_0.tsv.gz"), emit: input
    shell:
        '''
        concatenate_tsvs.py -o !{name}.tsv.gz !{tsvs}
        ln -s !{tsvs[0]} !{name}_in_0.tsv.gz # Link input to output for testing
        '''
}
