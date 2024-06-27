// Generate mapping between genome IDs and taxids
process MAKE_HUMAN_VIRUS_ID_MAPPING {
    label "biopython"
    label "single"
    input:
        path(ncbi_metadata)
        path(genbank_genomes)
    output:
        path("genomeid-to-taxid.json")
    shell:
        '''
        #!/usr/bin/env python
        # Import packages
        import json
        import gzip
        from Bio.SeqIO.FastaIO import SimpleFastaParser
        # Generate mapping
        genome_to_taxid = {}
        with(open("!{ncbi_metadata}") as inf):
            cols = None
            for line in inf:
                bits = line.rstrip("\\n").split("\\t")
                if not cols:
                    cols = bits
                    file_idx = cols.index("local_filename")
                    taxid_idx = cols.index("taxid")
                    continue
                with gzip.open(bits[file_idx], "rt") as inf2:
                    for title, sequence in SimpleFastaParser(inf2):
                        genome_id, name = title.split(" ", 1)
                        taxid = int(bits[taxid_idx])
                        genome_to_taxid[genome_id] = taxid, name
        # Write output
        out_path = "genomeid-to-taxid.json"
        with open(out_path, "w") as outf:
            json.dump(genome_to_taxid, outf)
        '''
}
