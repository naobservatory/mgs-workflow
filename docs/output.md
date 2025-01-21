# Outputs

If the pipeline runs to completion, the following output files are expected. In the future, we will add more specific information about the outputs, including in-depth descriptions of the columns in the output files.

All pipeline output can be found in the `output` directory, which is broken into four subdirectories:

- `input`: Directory containing saved input information (useful for trying to reproduce someone else's results)
- `logging`: Log files containing meta-level information about the pipeline run itself.
- `intermediates`: Intermediate files produced by key stages in the run workflow, saved for nonstandard downstream analysis.
- `results`: Directory containing processed results files for standard downstream analysis.

## Run workflow

Main heading represents the folder name, and subheadings represent a description of the file's usage. If the file is not in the heading folder name, the relative path is given.

### `input`

- `adapters.fasta`: FASTA file of adapter sequences used for adapter screening.
- `run-params.json`: JSON file giving all the parameters passed to the pipeline.
- `index-params.json`: JSON file giving parameters used to generate index directory (`params.ref_dir`).
- `samplesheet.csv`: Copy of the samplesheet file used to configure the pipeline (specified by `params.sample_sheet`).

### `logging`

- `pipeline-version.txt`: Version of the pipeline used for the run.
- `time.txt`: Start time of the run.
- `trace.txt`: Tab delimited log of all the information for each task run in the pipeline including runtime, memory usage, exit status, etc. Can be used to create an execution timeline using the the script `bin/plot-timeline-script.R` after the pipeline has finished running. More information regarding the trace file format can be found [here](https://www.nextflow.io/docs/latest/reports.html#trace-file).
- `pipeline-version-index.txt`: Version of pipeline used to generate index directory (`params.ref_dir`).

### `intermediates`

- `reads/raw_viral`: Directory containing raw reads corresponding to those reads that survive initial viral screening with BBDuk.

### `results`

#### `qc`
- `total_reads_qc.tsv.gz`: Total number of reads processed at each stage of the preprocessing phase (`stage`)
- `qc_basic_stats.tsv.gz`: Summary statistics for each sample at each stage of the preprocessing phase (`stage`), including:
    - GC content (`percent GC`);
    - Average read length (`mean_seq_len`);
    - Number of read pairs (`n_read pairs`);
    - Approximate number of base pairs in reads (`n_bases_approx`);
    - Percent duplicates as measured by FASTQC (`percent_duplicates`);
    - Pass/fail scores for each test conducted by FASTQC.
- `qc_adapter_stats.tsv.gz`: Adapter statistics calculated by FASTQC for each sample and preprocessing stage, given as a percentage of reads containing adapter content (`pc_adapters`) at each position along the read (`position`) for each adapter detected (`adapter`) for each read in the read pair (`read_pair`).
- `qc_quality_base_stats.tsv.gz`: Per-base read-quality statistics calculated by FASTQC for each sample and preprocessing stage, given as the mean Phred score (`mean_phred_score`) at each position along the read (`position`) for each read in the read pair (`read_pair`).
- `qc_quality_sequence_stats.tsv.gz`: Per-sequence read-quality statistics calculated by FASTQC for each sample and preprocessing stage, given as the number of reads (`n_sequences`) with a given mean Phred score (`mean_phred_score`) for each read in the read pair (`read_pair`).

#### `viral identification`
- `virus_hits_db.tsv.gz`: TSV output by the viral identification phase, giving information about each read pair assigned to a host-infecting virus.
- `virus_clade_counts.tsv.gz`: Summary of the previous file giving the number of HV read pairs mapped to each viral taxon. Includes both read pairs mapped directly to that taxon (`n_reads_direct`) and to that taxon plus all descendent taxa (`n_reads_clade`).

#### `taxonomic identification`
- `kraken_reports_merged.tsv.gz`: Kraken output reports in TSV format, labeled by sample and ribosomal status.
- `bracken_reports_merged.tsv.gz`: Bracken output reports in TSV format, labeled by sample and ribosomal status.

#### `blast`
- `blast_hits_paired.tsv.gz`: Summarized BLASTN output for putative HV read pairs, giving, for each read pair and subject taxid:
    - The number of reads in the read pair with high-scoring matches to that taxid (`n_reads`).
    - The best bitscores of alignments to that taxid for each matching read (`bitscore_max` and `bitscore_min`)[^bitscore].
- `blast_forward_only.tsv.gz`: Summarized BLASTN output for putative HV read pairs, giving, for each read pair and subject taxid:
- `blast_reverse_only.tsv.gz`: Summarized BLASTN output for putative HV read pairs, giving, for each read pair and subject taxid:

[^bitscore]: If only one read aligns to the target, these two fields will be identical. If not, they will give the higher and lower of the best bitscores for the two reads in the pair..


## Index workflow

Main heading represents the folder name, and subheadings describes the tool that consumes the file. Files that are consumed by multiple tools or are not consumed by any tools are put in the `General` subheading. If the file is not in the heading folder name, the relative path is given.

### `input`

- `index-params.json`: JSON file giving all the parameters passed to the pipeline (useful for trying to reproduce someone else's results).

### `logging`

- `pipeline-version.txt`: Version of the pipeline with which index directory was created.
- `time.txt`: Start time of index workflow run.
- `trace.txt`: Nextflow trace file containing logging information for each process performed during the workflow run.

### `results`

#### `General`

- `total-virus-db-annotated.tsv.gz`: Database generated from NCBI taxonomy and Virus-Host-DB giving taxonomy and host-infection information for each viral taxon.
- `taxonomy-nodes.dmp`: Taxonomy dump file from NCBI mapping between taxids and their parents in the NCBI taxonomy tree structure.
- `taxonomy-names.dmp`: Taxonomy dump file from NCBI mapping between taxids and taxon names.

#### `BLAST`

- `core_nt`: Directory containing extracted BLAST database files for BLASTing against core_nt.

#### `Bowtie2`

- `bt2-virus-index`: Directory containing Bowtie2 index for host-infecting viral genomes.
- `bt2-human-index`: Directory containing Bowtie2 index for the human genome.
- `bt2-other-index`: Directory containing Bowtie2 index for other contaminant sequences.
- `virus-genome-metadata-gid.tsv.gz`: Genome metadata file generated during download of HV genomes from viral Genbank, annotated additionally with Genome IDs used by Bowtie2 (allowing mapping between genome ID and taxid).

#### `Kraken2`

- `kraken_db`: Directory containing Kraken2 reference database (default: Most recent version of Standard).

#### `BBduk`

- `virus-genomes-filtered.fasta.gz`: FASTA file containing host-infecting viral genomes downloaded from viral Genbank (filtered to remove transgenic, contaminated, or erroneous sequences).
- `ribo-ref-concat.fasta.gz`: Reference database of ribosomal LSU and SSU sequences from SILVA.
