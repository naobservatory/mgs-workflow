// Collect full list of descendent taxids for each viral taxid
process GET_HV_DESCENDENTS {
    label "biopython"
    label "single"
    input:
        path(human_virus_db)
    output:
        path("human-virus-taxids-all.txt"), emit: taxids
        path("human-virus-taxid-descendents.json"), emit: descendents
    shell:
        '''
        #!/usr/bin/env python
        # Import packages
        from collections import defaultdict
        import subprocess
        import json
        # Load human viruses
        print("Importing virus db...", end="")
        human_viruses = {}
        with open("!{human_virus_db}") as inf:
            for line in inf:
                taxid, name = line.strip().split("\\t")
                human_viruses[int(taxid)] = name
        print("done.")
        # Get descendents of each viral taxid
        print("Fetching viral descendents:")
        taxid_descendents = defaultdict(list)
        taxid_unique = set()
        for hv_taxid in human_viruses:
            print("\tFetching descendants of", hv_taxid, human_viruses[hv_taxid], end="")
            taxid_unique.add(hv_taxid)
            taxid_descendents[hv_taxid].append(hv_taxid)
            # Get descendent taxids from NCBI
            cmd = ["gimme_taxa.py", str(hv_taxid)]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            if p.returncode == 1 and "taxid not found" in error.decode("utf-8"):
                print("Taxid not found: {}".format(hv_taxid))
            elif p.returncode != 0:
                raise Exception(error.decode("utf-8"))
            desc_split = output.decode("utf-8").split("\\n")
            # Add descendent taxids to lists
            tab = False # Skip frontmatter from stdout
            n_desc = 0
            for line in desc_split:
                line = line.strip()
                if not line:
                    continue
                elif line.startswith("parent_taxid"):
                    tab = True # Stop skipping after header line
                    continue
                elif not tab:
                    continue
                try:
                    parent_taxid, descendent_taxid, descendent_name = line.split("\\t")
                except ValueError:
                    print(line)
                    raise
                descendent_taxid = int(descendent_taxid)
                taxid_unique.add(descendent_taxid)
                taxid_descendents[hv_taxid].append(descendent_taxid)
                n_desc += 1
            print("-", n_desc, "total descendents fetched.")
        print("Fetching complete.")
        # Write output
        out_path_list = "human-virus-taxids-all.txt"
        out_path_json = "human-virus-taxid-descendents.json"
        with open(out_path_list, "w") as outf:
            for taxid in sorted(taxid_unique):
                outf.write("%s\\n" % taxid)
        with open(out_path_json, "w") as outf:
            json.dump(taxid_descendents, outf)
        '''
}
