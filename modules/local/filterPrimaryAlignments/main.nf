// Filter joined TSV to keep only primary alignments (is_secondary=False)
process FILTER_PRIMARY_ALIGNMENTS {
    label "python"
    label "single"
    
    input:
        tuple val(sample), path(joined_tsv)
    
    output:
        tuple val(sample), path("${sample}_primary_alignments.tsv.gz"), emit: filtered
    
    shell:
        '''
        set -o pipefail
        out=!{sample}_primary_alignments.tsv.gz
        zcat !{joined_tsv} | filter_primary_alignments.py | gzip > ${out}
        '''
}
