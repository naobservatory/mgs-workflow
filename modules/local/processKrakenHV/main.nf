// Process Kraken2 output and identify HV- and non-HV-assigned reads
process PROCESS_KRAKEN_HV {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(output)
        path nodes_path
        path hv_db_path
    output:
        tuple val(sample), path("${sample}_kraken_hv_processed.tsv")
    shell:
        '''
        in=!{output}
        out=!{sample}_kraken_hv_processed.tsv
        nodes=!{nodes_path}
        hv=!{hv_db_path}
        script_path=./!{script_process_kraken}
        process_kraken_hv.py ${in} ${hv} ${nodes} ${out}
        '''
}
