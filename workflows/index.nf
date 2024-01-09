/**************
| 0. PREAMBLE |
**************/

// Set up directories
pubDir = "${params.pub_dir}"

// Generate paths from param files
virus_db_path = "${params.virus_db}"

/*******************************************
| 1. RIBOSOMAL REFERENCE FOR RIBODEPLETION |
*******************************************/

// 1.1. Download & concatenate ribosomal references
process JOIN_RIBO_REF {
    label "BBTools"
    label "single"
    publishDir "${pubDir}/ribo", mode: "symlink"
    publishDir "${pubDir}/output", mode: "copy"
    input:
        val(ssu_url)
        val(lsu_url)
    output:
        path("ribo_ref_concat.fasta.gz")
    shell:
        '''
        # Download references
        wget !{ssu_url} -O ssu_ref.fasta.gz
        wget !{lsu_url} -O lsu_ref.fasta.gz
        in="ssu_ref.fasta.gz lsu_ref.fasta.gz"
        cat ${in} > ribo_ref_concat.fasta.gz # TODO: Make robust to differences in gzip status
        '''
}

workflow PREPARE_RIBOSOMAL {
    take:
        ssu_url
        lsu_url
    main:
        ref_ch = JOIN_RIBO_REF(ssu_url, lsu_url)
    emit:
        ref = ref_ch
}

/****************************************
| 2. HUMAN REFERENCE FOR HOST DEPLETION |
****************************************/

// 2.1. Obtain human genome reference
process GET_HUMAN {
    label "BBTools"
    label "single"
    publishDir "${pubDir}/human", mode: "symlink"
    input:
        val(human_url)
    output:
        path("human_ref.fasta.gz")
    shell:
        '''
        path=human_ref.fasta.gz
        wget !{human_url} -O ${path}
        '''
}

// 2.2. Mask human genome in preparation for indexing
process BBMASK_HUMAN {
    label "BBTools"
    label "large"
    publishDir "${pubDir}/human", mode: "symlink"
    input:
        path(reference_path)
    output:
        path("human_ref_masked.fasta.gz")
    shell:
        '''
        in=!{reference_path}
        out=human_ref_masked.fasta.gz
        io="in=${in} out=${out}"
        par="threads=!{task.cpus} maskrepeats masklowentropy"
        bbmask.sh ${par} ${io}
        '''
}

// 2.3. Build and tarball human index from reference genome
process BBMAP_INDEX_HUMAN {
    label "BBTools"
    label "large"
    publishDir "${pubDir}/human", mode: "symlink"
    publishDir "${pubDir}/output", mode: "copy"
    input:
        path(masked_reference)
    output:
        path("human_ref_index.tar.gz")
    shell:
        '''
        mkdir human_ref_index
        cp !{masked_reference} human_ref_index/human_ref.fasta.gz
        cd human_ref_index
        bbmap.sh ref=human_ref.fasta.gz t=!{task.cpus} -Xmx30g
        cd ..
        tar -czf human_ref_index.tar.gz human_ref_index
        '''
}

workflow PREPARE_HUMAN {
    take:
        human_url
    main:
        get_ch = GET_HUMAN(human_url)
        mask_ch = BBMASK_HUMAN(get_ch)
        index_ch = BBMAP_INDEX_HUMAN(mask_ch)
    emit:
        index = index_ch
}

/******************************************
| 3. OTHER REFERENCES FOR DECONTAMINATION |
******************************************/

// 3.1. Join reference sequences together
// TODO: Replace hardcoded references with an extensible list
process JOIN_OTHER_REF {
    label "BBTools"
    label "single"
    publishDir "${pubDir}/other", mode: "symlink"
    input:
        val(cow_url)
        val(pig_url)
    output:
        path("other_ref_concat.fasta.gz")
    shell:
        '''
        # Download references
        wget !{cow_url} -O cow_ref.fasta.gz
        wget !{pig_url} -O pig_ref.fasta.gz
        in="cow_ref.fasta.gz pig_ref.fasta.gz"
        cat ${in} > other_ref_concat.fasta.gz # TODO: Make robust to differences in gzip status
        '''
}

// 3.2. Mask reference genomes in preparation for indexing
process BBMASK_REFERENCES {
    label "BBTools"
    label "large"
    publishDir "${pubDir}/other", mode: "symlink"
    input:
        path(concat_references)
    output:
        path("other_ref_masked.fasta.gz")
    shell:
        '''
        in=!{concat_references}
        out=other_ref_masked.fasta.gz
        io="in=${in} out=${out}"
        par="threads=!{task.cpus} maskrepeats masklowentropy"
        bbmask.sh ${io} ${par}
        '''
}

// 3.3. Build index from reference genomes
process BBMAP_INDEX_REFERENCES {
    label "BBTools"
    label "max"
    publishDir "${pubDir}/other", mode: "symlink"
    publishDir "${pubDir}/output", mode: "copy"
    input:
        path(masked_reference)
    output:
        path("other_ref_index.tar.gz")
    shell:
        '''
        mkdir other_ref_index
        cp !{masked_reference} other_ref_index/other_ref.fasta.gz
        cd other_ref_index
        bbmap.sh ref=other_ref.fasta.gz t=!{task.cpus} usemodulo -Xmx60g
        cd ..
        tar -czf other_ref_index.tar.gz other_ref_index
        '''
}

workflow PREPARE_OTHER {
    take:
        cow_url
        pig_url
    main:
        join_ch = JOIN_OTHER_REF(cow_url, pig_url) // TODO: Replace hardcoded references with extensible list
        mask_ch = BBMASK_REFERENCES(join_ch)
        index_ch = BBMAP_INDEX_REFERENCES(mask_ch)
    emit:
        index = index_ch
}

// TODO: Add custom kraken DB construction?

/****************************************
| 4. HUMAN-VIRUS GENOMES FOR EXTRACTION |
****************************************/

// 4.1. Collect full list of descendent taxids for each viral taxid
process GET_HV_DESCENDENTS {
    label "Python"
    label "single"
    publishDir "${pubDir}/hviral", mode: "symlink"
    input:
        path(human_virus_path)
    output:
        path("human-virus-taxids-all.txt")
        path("human-virus-taxid-descendents.json")
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
        with open("!{human_virus_path}") as inf:
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

// 4.2. Get viral genome sequences based on expanded list of taxids
process GET_HV_GENOMES {
    label "Python"
    label "single"
    publishDir "${pubDir}/hviral", mode: "symlink"
    input:
        path(hv_taxids_all)
    output:
        path("ncbi_fetch_metadata.txt")
        path("genbank_genomes")
    shell:
        '''
        io="--taxids !{hv_taxids_all} --metadata-table ncbi_fetch_metadata.txt -o genbank_genomes"
        par="--section genbank --formats fasta --flat-output"
        ncbi-genome-download ${par} ${io} viral
        '''
}

// 4.3. Generate mapping between genome IDs and taxids
process MAKE_ID_MAPPING {
    label "Python"
    label "single"
    publishDir "${pubDir}/hviral", mode: "symlink"
    publishDir "${pubDir}/output", mode: "copy"
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

// 4.4. Collate HV genomes into a single file
process COLLATE_HV_GENOMES {
    label "BBTools"
    label "single"
    publishDir "${pubDir}/hviral", mode: "symlink"
    publishDir "${pubDir}/output", mode: "copy"
    input:
        path(genbank_genomes)
    output:
        path("human-viral-genomes.fasta.gz")
    shell:
        '''
        cat !{genbank_genomes}/*.fna.gz > human-viral-genomes.fasta.gz
        '''
}

// 4.5. Mask collated virus genomes
// TODO: Replace with bbmask?
process MASK_HV_GENOMES {
    label "BLAST"
    label "single"
    publishDir "${pubDir}/hviral/index", mode: "symlink"
    input:
        path(collected_genomes)
    output:
        path("masked_genomes.fasta.gz")
    shell:
        '''
        zcat -f !{collected_genomes} | dustmasker -out "masked_genomes.fasta" -outfmt fasta
        sed -i '/^>/!s/[a-z]/x/g' masked_genomes.fasta
        gzip masked_genomes.fasta
        '''
}

// 4.6. Build Bowtie2 index from masked genomes
process BUILD_BOWTIE2_DB {
    label "Bowtie2"
    label "max"
    publishDir "${pubDir}/hviral/index", mode: "symlink"
    publishDir "${pubDir}/output", mode: "copy"
    input:
        path(masked_genomes)
    output:
        path("bt2_hv_index.tar.gz") // Output directory for index files
    shell:
        '''
        mkdir bt2_hv_index
        bowtie2-build -f --threads !{task.cpus} !{masked_genomes} bt2_hv_index/hv_index
        tar -czf bt2_hv_index.tar.gz bt2_hv_index
        '''
}

workflow PREPARE_VIRAL {
    take:
        hv_path
    main:
        desc_ch = GET_HV_DESCENDENTS(hv_path)
        genome_ch = GET_HV_GENOMES(desc_ch[0])
        map_ch = MAKE_ID_MAPPING(genome_ch[0], genome_ch[1])
        collate_ch = COLLATE_HV_GENOMES(genome_ch[1])
        mask_ch = MASK_HV_GENOMES(collate_ch)
        index_ch = BUILD_BOWTIE2_DB(mask_ch)
    emit:
        index = index_ch
        map = map_ch
        genomes = collate_ch
}

/*********************************************************
| 5. PREPARE TAXONOMY REFERENCE FOR HUMAN VIRUS ANALYSIS |
*********************************************************/

// 5.1. Obtain taxonomy files from NCBI
process GET_TAXONOMY {
    label "BBTools"
    label "single"
    publishDir "${pubDir}/taxonomy", mode: "symlink"
    input:
        val(taxonomy_url)
    output:
        path("taxonomy.zip")
    shell:
        '''
        path=taxonomy.zip
        wget !{taxonomy_url} -O ${path}
        '''
}

// 5.2. Extract taxonomy archive and access nodes file
process EXTRACT_TAXONOMY {
    label "base"
    label "single"
    publishDir "${pubDir}/taxonomy", mode: "symlink"
    input:
        path(taxonomy_zip)
    output:
        path("taxonomy")
        path("taxonomy-nodes.dmp")
    shell:
        '''
        unzip !{taxonomy_zip} -d taxonomy
        cp taxonomy/nodes.dmp taxonomy-nodes.dmp
        '''
}

workflow PREPARE_TAXONOMY {
    take:
        tax_url
    main:
        get_ch = GET_TAXONOMY(tax_url)
        extract_ch = EXTRACT_TAXONOMY(get_ch)
    emit:
        dir = extract_ch[0]
        nodes = extract_ch[1]
}

/****************
| MAIN WORKFLOW |
****************/

workflow {
    PREPARE_RIBOSOMAL(params.ssu_url, params.lsu_url)
    PREPARE_HUMAN(params.human_url)
    PREPARE_OTHER(params.cow_url, params.pig_url)
    PREPARE_VIRAL(virus_db_path)
    PREPARE_TAXONOMY(params.taxonomy_url)
}
