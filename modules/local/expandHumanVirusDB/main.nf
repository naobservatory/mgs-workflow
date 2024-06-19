process EXPAND_HUMAN_VIRUS_DB {
    label "biopython"
    label "single"
    input:
        path(human_virus_db)
        path(taxonomy_nodes)
        path(taxonomy_names)
    output:
        path("human-virus-db-expanded.tsv")
    shell:
'''
#!/usr/bin/env python
import sys
# 1. Get raw HV taxids
raw_hv = set()
with open("!{human_virus_db}") as inf:
    for line in inf:
        taxid, name = line.strip().split("\\t")
        raw_hv.add(int(taxid))
# 2. Get parent/child taxid mapping from nodes
children = {}
with open("!{taxonomy_nodes}") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = line.replace("\\t|\\n", "").split("\\t|\\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        if parent_taxid not in children:
            children[parent_taxid] = []
        children[parent_taxid].append(child_taxid)
# 3. Get expanded list of HV taxids
hv = set()
def add_children(taxid):
    hv.add(taxid)
    for child in children.get(taxid, []):
        add_children(child)
for taxid in raw_hv:
    add_children(taxid)
# 4. Get taxonomic names
taxonomic_names = {}
with open("!{taxonomy_names}") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace("\\t|\\n", "").split("\\t|\\t")
        taxid = int(taxid)
        if taxid in hv:
            if taxid not in taxonomic_names or name_class == "scientific name":
                taxonomic_names[taxid] = name
# 5. Write output
with open("human-virus-db-expanded.tsv", "w") as outf:
    for taxid in sorted(hv):
        outf.write("%s\\t%s\\n" % (taxid, taxonomic_names[taxid]))
'''
}
