// Extract information from a paired Bowtie2 SAM alignment file into a TSV
process PROCESS_BOWTIE2_SAM_PAIRED {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(sam)
        path genomeid_taxid_map_path
    output:
        tuple val(sample), path("${sample}_bowtie2_sam_processed.tsv")
    shell:
        '''
        in=!{sam}
        out=!{sample}_bowtie2_sam_processed.tsv
        map=!{genomeid_taxid_map_path}
        process_bowtie2_sam.py ${in} ${map} ${out}
        '''
}
