// Collect full list of descendent taxids for each viral taxid
process GET_HV_DESCENDENTS {
    label "biopython"
    label "max"
    input:
        path(human_virus_db)
    output:
        path("human-virus-taxids-all.txt"), emit: taxids
        path("human-virus-taxid-descendents.json"), emit: descendents
    shell:
        '''
        get-hv-descendents.py !{human_virus_db} human-virus-taxids-all.txt human-virus-taxid-descendents.json --threads !{task.cpus}
        '''
}
