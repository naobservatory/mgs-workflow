// Label Kraken2 reports with sample IDs
process LABEL_KRAKEN_REPORTS {
    label "tidyverse"
    label "single"
    input:
        tuple val(sample), path(report)
    output:
        path("${sample}_labeled.report.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        cnames <- c("pc_reads_total", "n_reads_clade", "n_reads_direct",
                    "n_minimizers_total", "n_minimizers_distinct", "rank", "taxid", "name")
        tab <- read_tsv("!{report}", col_names = cnames, show_col_types = FALSE) %>%
            mutate(sample="!{sample}")
        write_tsv(tab, "!{sample}_labeled.report.tsv.gz")
        '''
}

// Label Kraken2 reports with sample IDs (streamed Python version)
process LABEL_KRAKEN_REPORTS_STREAMED {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(report)
    output:
        path("${sample}_labeled.report.tsv.gz"), emit: report
        path("${sample}_in.report.tsv.gz"), emit: input
    shell:
        '''
        # Run labeling script
        label_kraken_reports.py !{report} !{sample} !{sample}_labeled.report.tsv.gz
        # Link input files for testing
        ln -s !{report} !{sample}_in.report.tsv.gz
        '''
}
