pubDir = "${params.pub_dir}"

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { RAW } from "../subworkflows/local/raw" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "raw_concat")
include { CLEAN } from "../subworkflows/local/clean" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "cleaned")
include { DEDUP } from "../subworkflows/local/dedup" addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "dedup")
include { RIBODEPLETION as RIBO_INITIAL } from "../subworkflows/local/ribodepletion", addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribo_initial", min_kmer_fraction: "0.6", k: "43")
include { RIBODEPLETION as RIBO_SECONDARY } from "../subworkflows/local/ribodepletion", addParams(fastqc_cpus: "2", fastqc_mem: "4 GB", stage_label: "ribo_secondary", min_kmer_fraction: "0.4", k: "27")
include { TAXONOMY } from "../subworkflows/local/taxonomy"

/************************************
| 0. PREPARE REFERENCES AND INDEXES |
************************************/

// 0.1. Extract human index for Bowtie2
process EXTRACT_HUMAN_BOWTIE2 {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(human_index_tarball)
    output:
        path("bt2_human_index")
    shell:
        '''
        tar -xzf !{human_index_tarball}
        '''
}

// 0.2. Extract contaminant index for Bowtie2
process EXTRACT_OTHER_BOWTIE2 {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(other_index_tarball)
    output:
        path("bt2_other_index")
    shell:
        '''
        tar -xzf !{other_index_tarball}
        '''
}

// 0.1. Extract human index for BBMap
process EXTRACT_HUMAN_BBMAP {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(human_index_tarball)
    output:
        path("human_ref_index")
    shell:
        '''
        tar -xzf !{human_index_tarball}
        '''
}

// 0.2. Extract contaminant index for BBMap
process EXTRACT_OTHER_BBMAP {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(other_index_tarball)
    output:
        path("other_ref_index")
    shell:
        '''
        tar -xzf !{other_index_tarball}
        '''
}

// 0.3. Extract human-infecting virus index for Bowtie2
process EXTRACT_HV {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(hv_index_tarball)
    output:
        path("bt2_hv_index")
    shell:
        '''
        tar -xzf !{hv_index_tarball}
        '''
}

// 0.4. Extract Kraken DB
process EXTRACT_KRAKEN {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(kraken_tarball)
    output:
        path("kraken_db")
    shell:
        '''
        mkdir kraken_db
        tar -xzf !{kraken_tarball} -C kraken_db
        '''
}

// 0.5. Copy ribosomal reference
process COPY_RIBO {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(ribo_db)
    output:
        path("ribo_db_ref.fasta.gz")
    shell:
        '''
        cp !{ribo_db} ribo_db_ref.fasta.gz
        '''
}

// 0.6. Copy adapters
process COPY_ADAPTERS {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(adapters)
    output:
        path("adapters_ref.fasta")
    shell:
        '''
        cp !{adapters} adapters_ref.fasta
        '''
}

// 0.7. Copy viral taxa
process COPY_TAXA {
    label "BBTools"
    label "single"
    errorStrategy "retry"
    input:
        path(viral_taxa)
    output:
        path("viral-taxa.tsv.gz")
    shell:
        '''
        cp !{viral_taxa} viral-taxa.tsv.gz
        '''
}

// 0.8. Copy sample metadata CSV
process COPY_SAMPLE_METADATA {
    label "BBTools"
    label "single"
    publishDir "${pubDir}", mode: "copy", overwrite: "true"
    errorStrategy "retry"
    input:
        path(sample_metadata)
    output:
        path("sample-metadata.csv")
    shell:
        '''
        cp !{sample_metadata} sample-metadata.csv
        '''
}

workflow PREPARE_REFERENCES {
    take:
        human_index_bt2_path
        other_index_bt2_path
        human_index_bb_path
        other_index_bb_path
        hv_index_path
        kraken_db_path
        ribo_db_path
        adapter_path
        taxid_path
        sample_metadata_path
    main:
        human_ch_bt2 = EXTRACT_HUMAN_BOWTIE2(human_index_bt2_path)
        other_ch_bt2 = EXTRACT_OTHER_BOWTIE2(other_index_bt2_path)
        human_ch_bb = EXTRACT_HUMAN_BBMAP(human_index_bb_path)
        other_ch_bb = EXTRACT_OTHER_BBMAP(other_index_bb_path)
        hv_ch = EXTRACT_HV(hv_index_path)
        kraken_ch = EXTRACT_KRAKEN(kraken_db_path)
        ribo_ch = COPY_RIBO(ribo_db_path)
        adapter_ch = COPY_ADAPTERS(adapter_path)
        taxa_ch = COPY_TAXA(taxid_path)
        sample_ch = COPY_SAMPLE_METADATA(sample_metadata_path)
    emit:
        human_bt2 = human_ch_bt2
        other_bt2 = other_ch_bt2
        human_bb = human_ch_bb
        other_bb = other_ch_bb
        hv = hv_ch
        kraken = kraken_ch
        ribo = ribo_ch
        adapters = adapter_ch
        taxa = taxa_ch
        metadata = sample_ch
}

/****************************************
| 5. HUMAN VIRUS DETECTION WITH BOWTIE2 |
****************************************/

// 5.1. Run Bowtie2 and return mapped HV reads
process RUN_BOWTIE2 {
    label "Bowtie2"
    label "large"
    input:
        tuple val(sample), path(reads_noribo), path(reads_ribo), path(stats)
        path(index_dir)
    output:
        tuple val(sample), path("${sample}_bowtie2_mapped.sam"), path("${sample}_bowtie2_conc_{1,2}.fastq.gz"), path("${sample}_bowtie2_unconc_{1,2}.fastq.gz")
    shell:
        '''
        in1=!{reads_noribo[0]}
        in2=!{reads_noribo[1]}
        idx="!{index_dir}/hv_index"
        sam="!{sample}_bowtie2_mapped.sam"
        alc="!{sample}_bowtie2_conc_%.fastq.gz"
        unc="!{sample}_bowtie2_unconc_%.fastq.gz"
        io="-1 ${in1} -2 ${in2} -x ${idx} -S ${sam} --al-conc-gz ${alc} --un-conc-gz ${unc}"
        par="--threads !{task.cpus} --no-unal --no-sq --local --very-sensitive-local --score-min G,5,11"
        bowtie2 ${par} ${io}
        '''
}

// 5.2. Process Bowtie2 SAM output
// NB: Currently paired, need to update if switch to merged
process PROCESS_BOWTIE_SAM {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(sam_out), path(reads_conc), path(reads_unconc)
        path script_process_sam
        path genomeid_taxid_map_path
    output:
        tuple val(sample), path("${sample}_bowtie2_sam_processed.tsv")
    shell:
        '''
        in=!{sam_out}
        out=!{sample}_bowtie2_sam_processed.tsv
        map=!{genomeid_taxid_map_path}
        script_path=./!{script_process_sam}
        ${script_path} ${in} ${map} ${out}
        '''
}

// 5.3. Get list of singly or discordantly aligned read pairs (if any)
process GET_UNCONC_READ_IDS {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(sam_out), path(reads_conc), path(reads_unconc)
    output:
        tuple val(sample), path("${sample}_bowtie2_mapped_unconc.txt")
    shell:
        '''
        #!/usr/bin/env python
        import re
        # Define input and output paths
        inp = "!{sam_out}"
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
                    status_field = [f for f in line_split if "YT" in f][0]
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

// 5.4. Extract singly or discordantly aligned reads from Bowtie2 output
process EXTRACT_UNCONC_READS {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(sam_out), path(reads_conc), path(reads_unconc), path(id_file)
    output:
        tuple val(sample), path("${sample}_bowtie2_mapped_unconc_{1,2}.fastq.gz")
    shell:
        '''
        #!/usr/bin/env python
        from Bio import SeqIO
        import gzip
        # Define input and output paths
        inp0="!{id_file}"
        inp1="!{reads_unconc[0]}"
        inp2="!{reads_unconc[1]}"
        outp1="!{sample}_bowtie2_mapped_unconc_1.fastq.gz"
        outp2="!{sample}_bowtie2_mapped_unconc_2.fastq.gz"
        # Get list of ids
        print("Importing IDs....")
        with open(inp0, "r") as inf:
            ids = [line.strip() for line in inf.readlines()]
        print("Done. {} IDs imported.".format(len(ids)))
        # Filter forward reads
        print("Filtering forward reads...")
        with gzip.open(inp1, "rt") as inf, gzip.open(outp1, "wt") as outf:
            seqs = SeqIO.parse(inf, "fastq")
            for s in seqs:
                if s.id in ids:
                    SeqIO.write(s, outf, "fastq")
        # Filter reverse reads
        print("Filtering reverse reads...")
        with gzip.open(inp2, "rt") as inf, gzip.open(outp2, "wt") as outf:
            seqs = SeqIO.parse(inf, "fastq")
            for s in seqs:
                if s.id in ids:
                    SeqIO.write(s, outf, "fastq")
        '''
}

// 5.5. Combine concordantly and non-concordantly mapped read pairs
process COMBINE_MAPPED_READS {
    label "BBTools"
    label "single"
    input:
        tuple val(sample), path(sam_out), path(reads_conc), path(reads_unconc), path(reads_mapped_unconc)
    output:
        tuple val(sample), path("${sample}_bowtie2_mapped_all_{1,2}.fastq.gz")
    shell:
        '''
        inc1=!{reads_conc[0]}
        inc2=!{reads_conc[1]}
        inu1=!{reads_mapped_unconc[0]}
        inu2=!{reads_mapped_unconc[1]}
        out1=!{sample}_bowtie2_mapped_all_1.fastq.gz
        out2=!{sample}_bowtie2_mapped_all_2.fastq.gz
        cat $inc1 $inu1 > $out1
        cat $inc2 $inu2 > $out2
        '''
}

// 5.6. Segregate human reads with Bowtie2
process BOWTIE2_HUMAN_DEPLETION {
    label "large"
    label "Bowtie2"
    input:
        tuple val(sample), path(reads_mapped)
        path(index_dir)
    output:
        tuple val(sample), path("${sample}_bowtie2_human.sam"), path("${sample}_bowtie2_human_{1,2}.fastq.gz"), path("${sample}_bowtie2_nohuman_{1,2}.fastq.gz")
    shell:
        '''
        in1=!{reads_mapped[0]}
        in2=!{reads_mapped[1]}
        idx="!{index_dir}/human_index"
        sam="!{sample}_bowtie2_human.sam"
        alc="!{sample}_bowtie2_human_%.fastq.gz"
        unc="!{sample}_bowtie2_nohuman_%.fastq.gz"
        io="-1 ${in1} -2 ${in2} -x ${idx} -S ${sam} --al-conc-gz ${alc} --un-conc-gz ${unc}"
        par="--threads !{task.cpus} --local --very-sensitive-local"
        bowtie2 ${par} ${io}
        '''
}

process BBMAP_HUMAN_DEPLETION {
    label "large"
    label "BBTools"
    errorStrategy "retry"
    input:
        tuple val(sample), path(sam), path(reads_human), path(reads_nohuman)
        path(index_ref_dir)
    output:
        tuple val(sample), path("${sample}_bbmap_nohuman_{1,2}.fastq.gz"), path("${sample}_bbmap_human_{1,2}.fastq.gz"), path("${sample}_bbmap_human.stats.txt")
    shell:
        '''
        # Define input/output
        unmerged1=!{reads_nohuman[0]}
        unmerged2=!{reads_nohuman[1]}
        op1=!{sample}_bbmap_nohuman_1.fastq.gz
        op2=!{sample}_bbmap_nohuman_2.fastq.gz
        of1=!{sample}_bbmap_human_1.fastq.gz
        of2=!{sample}_bbmap_human_2.fastq.gz
        stats=!{sample}_bbmap_human.stats.txt
        io_unmerged="in=${unmerged1} in2=${unmerged2} outu=${op1} outu2=${op2} outm=${of1} outm2=${of2} statsfile=${stats} path=!{index_ref_dir}"
        # Define parameters
        par="minid=0.8 maxindel=4 bwr=0.25 bw=25 quickmatch minhits=2 t=!{task.cpus} -Xmx30g"
        # Execute
        bbmap.sh ${io_unmerged} ${par}
        '''
}

// 5.7. Segregate reference reads with Bowtie2
process BOWTIE2_REFERENCE_DEPLETION {
    label "large"
    label "Bowtie2"
    input:
        tuple val(sample), path(reads_nohuman), path(reads_human), path(stats)
        path(index_dir)
    output:
        tuple val(sample), path("${sample}_bowtie2_ref.sam"), path("${sample}_bowtie2_ref_{1,2}.fastq.gz"), path("${sample}_bowtie2_noref_{1,2}.fastq.gz")
    shell:
        '''
        in1=!{reads_nohuman[0]}
        in2=!{reads_nohuman[1]}
        idx="!{index_dir}/other_index"
        sam="!{sample}_bowtie2_ref.sam"
        alc="!{sample}_bowtie2_ref_%.fastq.gz"
        unc="!{sample}_bowtie2_noref_%.fastq.gz"
        io="-1 ${in1} -2 ${in2} -x ${idx} -S ${sam} --al-conc-gz ${alc} --un-conc-gz ${unc}"
        par="--threads !{task.cpus} --local --very-sensitive-local"
        bowtie2 ${par} ${io}
        '''
}

process BBMAP_REFERENCE_DEPLETION {
    label "BBTools"
    label "large"
    errorStrategy "retry"
    input:
        tuple val(sample), path(sam), path(reads_ref), path(reads_noref)
        path(index_ref_dir)
    output:
        tuple val(sample), path("${sample}_bbmap_noref_{1,2}.fastq.gz"), path("${sample}_bbmap_ref_{1,2}.fastq.gz"), path("${sample}_bbmap_other.stats.txt")
    shell:
        '''
        # Define input/output
        unmerged1=!{reads_noref[0]}
        unmerged2=!{reads_noref[1]}
        op1=!{sample}_bbmap_noref_1.fastq.gz
        op2=!{sample}_bbmap_noref_2.fastq.gz
        of1=!{sample}_bbmap_ref_1.fastq.gz
        of2=!{sample}_bbmap_ref_2.fastq.gz
        stats=!{sample}_bbmap_other.stats.txt
        io_unmerged="in=${unmerged1} in2=${unmerged2} outu=${op1} outu2=${op2} outm=${of1} outm2=${of2} statsfile=${stats} path=!{index_ref_dir}"
        # Define parameters
        par="minid=0.8 maxindel=4 bwr=0.25 bw=25 quickmatch minhits=2 t=!{task.cpus} usemodulo -Xmx30g"
        # Execute
        bbmap.sh ${io_unmerged} ${par}
        '''
}

// 5.10. Deduplicate merged HV candidate reads in an RC-sensitive manner
process DEDUP_HV {
    label "BBTools"
    label "large"
    errorStrategy "retry"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_bowtie2_mjc_dedup.fastq.gz")
    shell:
        '''
        # Define input/output
        in=!{reads}
        out=!{sample}_bowtie2_mjc_dedup.fastq.gz
        io="in=${in} out=${out}"
        # Define parameters
        par="reorder dedupe containment t=!{task.cpus} -Xmx15g"
        # Execute
        clumpify.sh ${io} ${par}
        '''
}

// 5.12. Process Kraken2 output and identify HV- and non-HV-assigned reads
process PROCESS_KRAKEN_BOWTIE {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(output), path(report)
        path script_process_kraken
        path nodes_path
        path hv_path
    output:
        tuple val(sample), path("${sample}_bowtie2_kraken_processed.tsv")
    shell:
        '''
        in=!{output}
        out=!{sample}_bowtie2_kraken_processed.tsv
        nodes=!{nodes_path}
        hv=!{hv_path}
        script_path=./!{script_process_kraken}
        ${script_path} ${in} ${hv} ${nodes} ${out}
        '''
}

// 5.13. Merge processed SAM and Kraken TSVs and compute length-normalized alignment scores
process MERGE_SAM_KRAKEN {
    label "tidyverse"
    label "single"
    input:
        tuple val(sample), path(kraken_processed), path(sam_processed)
    output:
        path("${sample}_hv_hits_putative.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        sam <- read_tsv("!{sam_processed}", show_col_types = FALSE)
        krk <- read_tsv("!{kraken_processed}", show_col_types = FALSE)
        mrg <- sam %>% rename(seq_id = query_name) %>% inner_join(krk, by="seq_id")
        cat(nrow(sam), nrow(krk), nrow(mrg), "\n")
        cat(ncol(sam), ncol(krk), ncol(mrg), "\n")
        cat(names(mrg), "\n")
        cat(head(mrg$best_alignment_score_fwd))
        cat(head(mrg$query_len_fwd))
        cat(class(mrg$query_len_fwd))
        if (nrow(mrg)>0) {
            mrg <- mrg %>%
                mutate(adj_score_fwd = best_alignment_score_fwd/log(query_len_fwd),
                       adj_score_rev = best_alignment_score_rev/log(query_len_rev),
                       sample="!{sample}")
        }
        write_tsv(mrg, "!{sample}_hv_hits_putative.tsv.gz")
        '''
}

// 5.15. Perform initial HV read filtering
process FILTER_HV {
    label "tidyverse"
    label "single"
    input:
        path(hv_hits)
        val(score_threshold)
    output:
        path("hv_hits_putative_filtered.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        score_threshold <- !{score_threshold}
        data <- read_tsv("!{hv_hits}", col_names = TRUE, show_col_types = FALSE)
        filtered <- mutate(data, hit_hv = as.logical(!is.na(str_match(encoded_hits, paste0(" ", as.character(taxid), ":"))))) %>%
            mutate(adj_score_fwd = replace_na(adj_score_fwd, 0), adj_score_rev = replace_na(adj_score_rev, 0)) %>%
            filter((!classified) | assigned_hv) %>% 
            filter(adj_score_fwd > score_threshold | adj_score_rev > score_threshold | assigned_hv)
        print(dim(data))
        print(dim(filtered))
        write_tsv(filtered, "hv_hits_putative_filtered.tsv.gz")
        '''
}

// Collapse separate read pair entries in HV DB
process COLLAPSE_HV {
    label "tidyverse"
    label "single"
    publishDir "${pubDir}", mode: "copy", overwrite: "true"
    input:
        path(hv_hits_filtered)
    output:
        path("hv_hits_putative_collapsed.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        rmax <- function(x){
            if (all(is.na(x))) return(NA)
            return(max(x, na.rm = TRUE))
        }
        collapse <- function(x) ifelse(all(x == x[1]), x[1], paste(x, collapse="/"))
        hits_filtered <- read_tsv("!{hv_hits_filtered}", col_names = TRUE, show_col_types = FALSE)
        print(dim(hits_filtered))
        reads_collapsed <- hits_filtered %>% group_by(seq_id) %>% summarize(
            sample = collapse(sample), genome_id = collapse(genome_id),
            taxid_best = taxid[1], taxid = collapse(as.character(taxid)),
            best_alignment_score_fwd = rmax(best_alignment_score_fwd),
            best_alignment_score_rev = rmax(best_alignment_score_rev),
            query_len_fwd = rmax(query_len_fwd), query_seq_fwd = query_seq_fwd[!is.na(query_seq_fwd)][1],
            query_len_rev = rmax(query_len_rev), query_seq_rev = query_seq_rev[!is.na(query_seq_rev)][1],
            classified = rmax(classified), assigned_name = collapse(assigned_name),
            assigned_taxid_best = assigned_taxid[1], assigned_taxid = collapse(as.character(assigned_taxid)),
            assigned_hv = rmax(assigned_hv), hit_hv = rmax(hit_hv), encoded_hits = collapse(encoded_hits),
            adj_score_fwd = rmax(adj_score_fwd), adj_score_rev = rmax(adj_score_rev)
            ) %>% mutate(adj_score_max = pmax(adj_score_fwd, adj_score_rev))
        print(dim(reads_collapsed))
        write_tsv(reads_collapsed, "hv_hits_putative_collapsed.tsv.gz")
        '''
}

// 5.16. Extract FASTA from filtered sequences
process MAKE_HV_FASTA {
    label "tidyverse"
    label "single"
    input:
        path(hv_hits_collapsed)
    output:
        path("hv_hits_putative_collapsed.fasta")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        tab <- read_tsv("!{hv_hits_collapsed}", col_names = TRUE, show_col_types = FALSE)
        fasta_tab <- mutate(tab, seq_head = paste0(">", seq_id),
                            header1 = paste0(seq_head, "_1"), header2 = paste0(seq_head, "_2")) %>%
            select(header1, seq1=query_seq_fwd, header2, seq2=query_seq_rev)
        fasta_out <- do.call(paste, c(fasta_tab, sep="\n")) %>% paste(collapse="\n")
        write(fasta_out, "hv_hits_putative_collapsed.fasta")
        '''
}

// 5.17. Extract table of clade counts from HV reads
process COUNT_HV_CLADES {
    label "R"
    label "large"
    publishDir "${pubDir}", mode: "copy", overwrite: "true"
    input:
        path(hv_hits_collapsed)
        path(viral_taxa)
        path(script_count_taxa)
    output:
        path("hv_clade_counts.tsv.gz")
    shell:
        '''
        script_path=./!{script_count_taxa}
        ${script_path} --reads !{hv_hits_collapsed} --taxa !{viral_taxa} --output hv_clade_counts.tsv.gz
        '''
}

workflow MAP_HUMAN_VIRUSES {
    take:
        ribo_ch
        bt2_index_ch
        kraken_db_ch
        nodes_path
        hv_db_path
        genomeid_map_path
        human_index_bt2_ch
        other_index_bt2_ch
        human_index_bb_ch
        other_index_bb_ch
        viral_taxa_ch
        script_process_sam
        script_join_fastq
        script_process_kraken
        script_count_taxa
        bt2_score_threshold
        adapters
    main:
        // Run Bowtie2 and process output
        bowtie2_ch = RUN_BOWTIE2(ribo_ch, bt2_index_ch)
        bowtie2_processed_ch = PROCESS_BOWTIE_SAM(bowtie2_ch, script_process_sam, genomeid_map_path)
        bowtie2_unconc_ids_ch = GET_UNCONC_READ_IDS(bowtie2_ch)
        bowtie2_unconc_reads_ch = EXTRACT_UNCONC_READS(bowtie2_ch.combine(bowtie2_unconc_ids_ch, by: 0))
        bowtie2_reads_combined_ch = COMBINE_MAPPED_READS(bowtie2_ch.combine(bowtie2_unconc_reads_ch, by: 0))
        // Filter contaminants
        human_bt2_ch = BOWTIE2_HUMAN_DEPLETION(bowtie2_reads_combined_ch, human_index_bt2_ch)
        human_bb_ch = BBMAP_HUMAN_DEPLETION(human_bt2_ch, human_index_bb_ch)
        other_bt2_ch = BOWTIE2_REFERENCE_DEPLETION(human_bb_ch, other_index_bt2_ch)
        other_bb_ch = BBMAP_REFERENCE_DEPLETION(other_bt2_ch, other_index_bb_ch)
        // Merge & join for Kraken input
        merge_ch = BBMERGE_BOWTIE(other_bb_ch)
        join_ch = JOIN_BOWTIE(merge_ch, script_join_fastq)
        // Deduplicate
        dedup_ch = DEDUP_HV(join_ch)
        // Run Kraken2 and process output
        kraken_ch = KRAKEN_BOWTIE(dedup_ch, kraken_db_ch)
        kraken_processed_ch = PROCESS_KRAKEN_BOWTIE(kraken_ch, script_process_kraken, nodes_path, hv_db_path)
        merged_input_ch = kraken_processed_ch.combine(bowtie2_processed_ch, by: 0)
        // Merge and filter output
        merged_ch = MERGE_SAM_KRAKEN(merged_input_ch)
        merged_ch_2 = MERGE_SAMPLES_HV(merged_ch.collect().ifEmpty([]))
        filtered_ch = FILTER_HV(merged_ch_2, bt2_score_threshold)
        collapsed_ch = COLLAPSE_HV(filtered_ch)
        fasta_ch = MAKE_HV_FASTA(collapsed_ch)
        // Count clades
        count_ch = COUNT_HV_CLADES(collapsed_ch, viral_taxa_ch, script_count_taxa)
    emit:
        data_all = merged_ch_2
        data_collapsed = collapsed_ch
        fasta = fasta_ch
        counts = count_ch
}

/***************************************
| 6. BLAST VALIDATION (OPTIONAL, SLOW) |
***************************************/

// 6.0. Subsequence merged reads for BLAST with seqtk
process SUBSET_BLAST {
    label "seqtk"
    label "single"
    input:
        path(hv_hits_collapsed_fasta)
        val readFraction
    output:
        path("hv_hits_putative_subset.fasta")
    shell:
        '''
        # Define input/output
        in1=!{hv_hits_collapsed_fasta}
        out1=hv_hits_putative_subset.fasta
        # Count reads for validation
        echo "Input reads: $(zcat ${in1} | wc -l | awk '{ print $1/2 }')"
        # Carry out subsetting
        seed=${RANDOM}
        seqtk sample -s ${seed} ${in1} !{readFraction} > ${out1}
        # Count reads for validation
        echo "Output reads: $(wc -l ${out1} | awk '{ print $1/2 }')"
        '''
}

// 6.1. BLAST putative HV hits against nt
process BLAST_HV {
    label "BLAST"
    cpus 32
    memory "256 GB"
    input:
        path(hv_hits_subset_fasta)
        path(blast_nt_dir)
    output:
        path("hv_hits_putative_subset.blast.gz")
    shell:
        '''
        # Specify parameters
        infile=!{hv_hits_subset_fasta}
        db=!{blast_nt_dir}/nt
        outfile=hv_hits_putative_subset.blast
        threads=!{task.cpus}
        # Set up command
        io="-query ${infile} -out ${outfile} -db ${db}"
        par="-perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30 -num_threads ${threads}"
        # Run BLAST
        blastn ${io} ${par} -outfmt "6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send"
        # Gzip output
        gzip hv_hits_putative_subset.blast
        '''
}

// 6.2. Process & filter BLAST output
process FILTER_BLAST {
    label "tidyverse"
    cpus 1
    memory "16 GB"
    input:
        path(blast_results)
    output:
        path("hv_hits_blast_filtered.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        blast_results_path <- "!{blast_results}"
        blast_cols <- c("qseqid", "sseqid", "sgi", "staxid", "qlen", "evalue",
            "bitscore", "qcovs", "length", "pident", "mismatch", "gapopen",
            "sstrand", "qstart", "qend", "sstart", "send")
        blast_results <- read_tsv(blast_results_path, show_col_types = FALSE,
            col_names = blast_cols, col_types = cols(.default="c"))
        # Filter for best hit for each query/subject combination
        blast_results_best <- blast_results %>% group_by(qseqid, staxid) %>% 
          filter(bitscore == max(bitscore)) %>%
          filter(length == max(length)) %>% filter(row_number() == 1)
        # Rank hits for each query and filter for high-ranking hits
        blast_results_ranked <- blast_results_best %>% 
          group_by(qseqid) %>% mutate(rank = dense_rank(desc(bitscore)))
        blast_results_highrank <- blast_results_ranked %>% filter(rank <= 5) %>%
            mutate(read_pair = str_split(qseqid, "_") %>% sapply(nth, n=-1),
                   seq_id = str_split(qseqid, "_") %>% sapply(nth, n=1))
        # Write output
        write_tsv(blast_results_highrank, "hv_hits_blast_filtered.tsv.gz")
        '''
}

// 6.3. Collapse results across read pairs
process PAIR_BLAST {
    label "tidyverse"
    label "single"
    publishDir "${pubDir}", mode: "copy", overwrite: "true"
    input:
        path(blast_results_filtered)
    output:
        path("hv_hits_blast_paired.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        blast_results_filtered_path <- "!{blast_results_filtered}"
        blast_results_filtered <- read_tsv(blast_results_filtered_path, show_col_types = FALSE)
        # Summarize by read pair and taxid
        blast_results_paired <- blast_results_filtered %>%
            mutate(bitscore = as.numeric(bitscore)) %>%
            group_by(seq_id, staxid) %>%
            summarize(bitscore_max = max(bitscore), bitscore_min = min(bitscore),
                      n_reads = n(), .groups = "drop")
        # Write output
        write_tsv(blast_results_paired, "hv_hits_blast_paired.tsv.gz")
        '''
}

workflow BLAST_HUMAN_VIRUSES {
    take:
        hv_hits_collapsed_fasta
        blast_nt_dir
        readFraction
    main:
        // Subset HV reads for BLAST
        subset_ch = SUBSET_BLAST(hv_hits_collapsed_fasta, readFraction)
        // BLAST putative HV hits against nt
        blast_ch = BLAST_HV(subset_ch, blast_nt_dir)
        // Process output
        filter_ch = FILTER_BLAST(blast_ch)
        pair_ch = PAIR_BLAST(filter_ch)
    emit:
        blast = blast_ch
        filtered = filter_ch
        paired = pair_ch
}

/**********************************
| 8. COLLATE AND PROCESS RESULTS |
**********************************/

// 8.1. Extract MultiQC data into more usable forms
process SUMMARIZE_MULTIQC_SINGLE {
    label "R"
    label "single"
    input:
        tuple val(stage), path(multiqc_data)
        path(script_summarize_multiqc)
    output:
        path("${stage}_qc_basic_stats.tsv.gz")
        path("${stage}_qc_adapter_stats.tsv.gz")
        path("${stage}_qc_quality_base_stats.tsv.gz")
        path("${stage}_qc_quality_sequence_stats.tsv.gz")
    shell:
        '''
        script_path=./!{script_summarize_multiqc}
        ${script_path} -i !{multiqc_data} -s !{stage} -o ${PWD}
        '''
}

// 8.2. Combine MultiQC summary files across workflow stages
process MERGE_MULTIQC {
    label "tidyverse"
    label "single"
    publishDir "${pubDir}", mode: "copy", overwrite: "true"
    input:
        path(basic_stats_tsvs)
        path(adapter_tsvs)
        path(base_quality_tsvs)
        path(sequence_quality_tsvs)
    output:
        path("qc_basic_stats.tsv.gz")
        path("qc_adapter_stats.tsv.gz")
        path("qc_quality_base_stats.tsv.gz")
        path("qc_quality_sequence_stats.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        # Import data
        in_paths_basic <- str_split("!{basic_stats_tsvs}", " ")[[1]]
        in_paths_adapt <- str_split("!{adapter_tsvs}", " ")[[1]]
        in_paths_qbase <- str_split("!{base_quality_tsvs}", " ")[[1]]
        in_paths_qseqs <- str_split("!{sequence_quality_tsvs}", " ")[[1]]
        tabs_basic <- lapply(in_paths_basic, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_adapt <- lapply(in_paths_adapt, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_qbase <- lapply(in_paths_qbase, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_qseqs <- lapply(in_paths_qseqs, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        # Bind rows
        tab_basic_out <- bind_rows(tabs_basic)
        tab_adapt_out <- bind_rows(tabs_adapt)
        tab_qbase_out <- bind_rows(tabs_qbase)
        tab_qseqs_out <- bind_rows(tabs_qseqs)
        # Write output
        write_tsv(tab_basic_out, "qc_basic_stats.tsv.gz")
        write_tsv(tab_adapt_out, "qc_adapter_stats.tsv.gz")
        write_tsv(tab_qbase_out, "qc_quality_base_stats.tsv.gz")
        write_tsv(tab_qseqs_out, "qc_quality_sequence_stats.tsv.gz")
        '''
}

// 8.3. Summarize taxonomic composition from Bracken and MultiQC output
process SUMMARIZE_COMPOSITION {
    label "tidyverse"
    label "single"
    publishDir "${pubDir}", mode: "copy", overwrite: "true"
    input:
        path(bracken_merged)
        path(multiqc_basic_merged)
    output:
        path("taxonomic_composition.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        # Import data
        bracken <- read_tsv("!{bracken_merged}", show_col_types = FALSE)
        basic   <- read_tsv("!{multiqc_basic_merged}", show_col_types = FALSE)
        total_assigned <- bracken %>% group_by(sample) %>% summarize(
          name = "Total",
          kraken_assigned_reads = sum(kraken_assigned_reads),
          added_reads = sum(added_reads),
          new_est_reads = sum(new_est_reads),
          fraction_total_reads = sum(fraction_total_reads)
        )
        bracken_spread <- bracken %>% select(name, sample, new_est_reads) %>%
          mutate(name = tolower(name)) %>%
          pivot_wider(id_cols = "sample", names_from = "name", values_from = "new_est_reads",
                      values_fill = 0)
        # Count reads
        read_counts_preproc <- basic %>% select(sample, stage, n_read_pairs) %>%
          pivot_wider(id_cols = c("sample"), names_from="stage", values_from="n_read_pairs",
                      values_fill = 0)
        read_counts <- read_counts_preproc %>%
          inner_join(total_assigned %>% select(sample, new_est_reads), by = "sample") %>%
          rename(assigned = new_est_reads) %>%
          inner_join(bracken_spread, by="sample")
        # Assess composition
        read_comp <- transmute(read_counts, sample=sample,
                               n_filtered = raw_concat-cleaned,
                               n_duplicate = cleaned-dedup,
                               n_ribosomal = (dedup-ribo_initial) + (ribo_initial-ribo_secondary),
                               n_unassigned = ribo_secondary-assigned,
                               n_bacterial = bacteria,
                               n_archaeal = archaea,
                               n_viral = viruses,
                               n_human = eukaryota)
        read_comp_long <- pivot_longer(read_comp, -(sample), names_to = "classification",
                                       names_prefix = "n_", values_to = "n_reads") %>%
          mutate(classification = fct_inorder(str_to_sentence(classification))) %>%
          group_by(sample) %>% mutate(p_reads = n_reads/sum(n_reads))
        # Write output
        write_tsv(read_comp_long, "taxonomic_composition.tsv.gz")
        '''
}

workflow PROCESS_OUTPUT {
    take:
        multiqc_data_raw_concat
        multiqc_data_cleaned
        multiqc_data_dedup
        multiqc_data_ribo_initial
        multiqc_data_ribo_secondary
        bracken_merged
        script_summarize_multiqc
    main:
        // Summarize each MultiQC directory separately
        multiqc_single = multiqc_data_raw_concat.mix(multiqc_data_cleaned, multiqc_data_dedup, multiqc_data_ribo_initial, multiqc_data_ribo_secondary)
        multiqc_summ   = SUMMARIZE_MULTIQC_SINGLE(multiqc_single, script_summarize_multiqc)
        multiqc_basic = multiqc_summ[0].collect().ifEmpty([])
        multiqc_adapt = multiqc_summ[1].collect().ifEmpty([])
        multiqc_qbase = multiqc_summ[2].collect().ifEmpty([])
        multiqc_qseqs = multiqc_summ[3].collect().ifEmpty([])
        // Merge MultiQC outputs
        multiqc_merged = MERGE_MULTIQC(multiqc_basic, multiqc_adapt, multiqc_qbase, multiqc_qseqs)
        // Summarize taxonomic composition
        taxo = SUMMARIZE_COMPOSITION(bracken_merged, multiqc_merged[0])
    emit:
        basic = multiqc_merged[0]
        adapt = multiqc_merged[1]
        qbase = multiqc_merged[2]
        qseqs = multiqc_merged[3]
        composition = taxo
}

/*****************
| MAIN WORKFLOWS |
*****************/

// Prepare libraries
libraries_ch = Channel
    .fromPath(params.library_tab)
    .splitCsv(header: true)
    .map{row -> [row.sample, row.library]}
    .groupTuple()

// Run first on new data to identify adapters
workflow prelim {
    PREPARE_REFERENCES(params.human_index_bt2, params.other_index_bt2, params.human_index_bb, params.other_index_bb, params.hv_index, params.kraken_db, params.ribo_db, params.adapters, params.viral_taxa_db, params.sample_tab)
    HANDLE_RAW_READS(libraries_ch, params.raw_dir)
    CLEAN_READS(HANDLE_RAW_READS.out.data, params.adapters)
}

// Complete primary workflow
workflow {
    // Prepare references & indexes
    PREPARE_REFERENCES(params.human_index_bt2, params.other_index_bt2, params.human_index_bb, params.other_index_bb, params.hv_index, params.kraken_db, params.ribo_db, params.adapters, params.viral_taxa_db, params.sample_tab)
    // Preprocessing
    RAW(libraries_ch, params.raw_dir, params.truncate_reads, params.n_reads_trunc)
    CLEAN(RAW.out.reads, PREPARE_REFERENCES.out.adapters)
    DEDUP(CLEAN_READS.out.reads)
    RIBO_INITIAL(DEDUP.out.reads, PREPARE_REFERENCES.out.ribo)
    RIBO_SECONDARY(RIBO_INITIAL.out.reads, PREPARE_REFERENCES.out.ribo)
    // Taxonomic profiling (all ribodepleted)
    TAX = TAXONOMY(RIBO_SECONDARY.out.reads, PREPARE_REFERENCES.out.kraken, 1, false, "D")
    // Taxonomic profiling (pre/post dedup)
    TAX_PRE = TAXONOMY(CLEAN.out.reads, PREPARE_REFERENCES.out.kraken, params.classify_dedup_subset, false, "D")
    TAX_POST = TAXONOMY(DEDUP.out.reads, PREPARE_REFERENCES.out.kraken, params.classify_dedup_subset, false, "D")
    
    // Human viral reads
    MAP_HUMAN_VIRUSES(REMOVE_RIBO_INITIAL.out.data, PREPARE_REFERENCES.out.hv, PREPARE_REFERENCES.out.kraken,
        params.nodes, params.hv_db, params.genomeid_map, 
        PREPARE_REFERENCES.out.human_bt2, PREPARE_REFERENCES.out.other_bt2, PREPARE_REFERENCES.out.human_bb, PREPARE_REFERENCES.out.other_bb,
        PREPARE_REFERENCES.out.taxa,
        params.script_process_sam, params.script_join_fastq, params.script_process_kraken, params.script_count_taxa,
        params.bt2_score_threshold, PREPARE_REFERENCES.out.adapters)
    if ( params.blast_hv ) {
        BLAST_HUMAN_VIRUSES(MAP_HUMAN_VIRUSES.out.fasta, params.blast_nt_dir, params.blast_fraction)
    }
    // Process output
    PROCESS_OUTPUT(HANDLE.multiqc_data, CLEAN_READS.out.multiqc_data, DEDUP_READS.out.multiqc_data, REMOVE_RIBO_INITIAL.out.multiqc_data, REMOVE_RIBO_SECONDARY.out.multiqc_data, CLASSIFY_READS.out.bracken, params.script_summarize_multiqc)
}
