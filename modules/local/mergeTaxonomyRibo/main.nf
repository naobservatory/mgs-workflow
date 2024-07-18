// Label and merge ribosomal and nonribosomal taxonomy profiling data
process MERGE_TAXONOMY_RIBO {
    label "tidyverse"
    label "single"
    input:
        path(kraken_reports_ribo)
        path(kraken_reports_noribo)
        path(bracken_ribo)
        path(bracken_noribo)
    output:
        path("kraken_reports_merged.tsv.gz"), emit: kraken
        path("bracken_reports_merged.tsv.gz"), emit: bracken
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        # Import data
        kr_ribo <- read_tsv("!{kraken_reports_ribo}", col_types = cols(.default="c"))
        kr_noribo <- read_tsv("!{kraken_reports_noribo}", col_types = cols(.default="c"))
        br_ribo <- read_tsv("!{bracken_reports_ribo}", col_types = cols(.default="c"))
        br_noribo <- read_tsv("!{bracken_reports_noribo}", col_types = cols(.default="c"))
        # Label and merge
        kr_out <- kr_ribo %>% mutate(ribosomal = TRUE) %>%
            bind_rows(kr_noribo %>% mutate(ribosomal = FALSE))
        br_out <- br_ribo %>% mutate(ribosomal = TRUE) %>%
            bind_rows(br_noribo %>% mutate(ribosomal = FALSE))
        # Write output
        write_tsv(kr_out, "kraken_reports_merged.tsv.gz")
        write_tsv(br_out, "bracken_reports_merged.tsv.gz")
        '''
}
