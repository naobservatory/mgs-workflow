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

// Label Bracken reports with sample IDs (streamed Python version)
process LABEL_BRACKEN_REPORTS_STREAMED {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(report)
    output:
        path("${sample}_labeled.bracken.tsv.gz"), emit: report
        path("${sample}_in.bracken.tsv.gz"), emit: input
    shell:
        '''
        # Run labeling script
        label_bracken_reports.py !{report} !{sample} !{sample}_labeled.bracken.tsv.gz
        # Link input files for testing
        ln -s !{report} !{sample}_in.bracken.tsv.gz
        '''
}
