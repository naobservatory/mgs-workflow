// From a Bowtie2 SAM file, get a list of IDs for singly or discordantly aligned read pairs (if any)
process EXTRACT_UNCONC_READ_IDS {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(sam)
    output:
        tuple val(sample), path("${sample}_bowtie2_mapped_unconc.txt")
    shell:
        '''
        #!/usr/bin/env python
        import re
        # Define input and output paths
        inp = "!{sam}"
        outp = "!{sample}_bowtie2_mapped_unconc.txt"
        # Extract read IDs from SAM file
        ids_out = set()
        with open(inp, "r") as inf:
            for line in inf:
                line_split = line.strip().split("\\t")
                if line_split[0].startswith("@"):
                    continue
                id = line_split[0]
                try:
                    status_field = [f for f in line_split if "YT:Z:" in f][0]
                    status = re.findall("YT:Z:(.*)", status_field)[0]
                except:
                    print(line_split)
                    raise
                if status in ["DP", "UP"]:
                    ids_out.add(id)
        # Write to new file
        ids_out_write = "\\n".join(list(ids_out)) + "\\n"
        with open(outp, "w") as outf:
            outf.write(ids_out_write)
        '''
}
