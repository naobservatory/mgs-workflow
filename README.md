# Nucleic Acid Observatory Viral Metagenomics Pipeline

This Nextflow pipeline is designed to process metagenomic sequencing data, characterize overall taxonomic composition, and identify and quantify reads mapping to viruses infecting certain host taxa of interest. It was developed as part of the [Nucleic Acid Observatory](https://naobservatory.org/) project.

## Pipeline description

### Overview

The pipeline currently consists of two workflows, an **index** workflow and a **run** workflow. The former creates index and reference files that are used by the latter, and is intended to be run first, after which many instantiations of the latter workflow can use the same index output files.

The run workflow then consists of five phases:

1. A **preprocessing phase**, in which input files undergo adapter & quality trimming.
2. A **viral identification phase**, in which a custom multi-step pipeline based around BBDuk, Bowtie2 and Kraken2 is used to sensitively and specifically identify reads from viruses infecting specific host taxa of interest ("host-infecting virus reads", or HV reads).
3. A **taxonomic profiling phase**, in which a random subset of reads (default 1M/sample) undergo ribosomal classification with BBDuk, followed by broader taxonomic classification with Kraken2.
4. An optional **BLAST validation phase**, in which putative HV reads from phase 2 are checked against nt to evaluate the sensitivity and specificity of the HV identification process.
5. A final **QC and output phase**, in which FASTQC, MultiQC and other tools are used to assess the quality of the data produced by the pipeline at various steps and produce summarized output data for downstream analysis and visualization.

A slight simplification of the overall process is given by this diagram:

![A flowchart summarizing the pipeline's run workflow.](/readme-workflow-diagram.png)

### Index workflow

The index workflow generates static index files from reference data. These reference data and indices don't depend on the reads processed by the run workflow, so a set of indices built by the index workflow can be used by multiple instantiations of the run workflow — no need to rebuild indices for every set of reads. The index workflow should be run once per reasonable length of time, balancing two considerations: Many of these index/reference files are derived from publicly available reference genomes or other resources, and should thus be updated and re-run periodically as new versions of these become available; however, to keep results comparable across datasets analyzed with the run workflow, this should be done relatively rarely.

Key inputs to the index workflow include:
- A TSV specifying names and taxonomic IDs (taxids) for host taxa for which to search for host-infecting viruses.
- A URL linking to a suitable Kraken2 database for taxonomic profiling (typically the [latest release](https://benlangmead.github.io/aws-indexes/k2) of the `k2_standard` database).
- URLS for up-to-date releases of reference genomes for various common contaminant species that can confound the identification of HV reads (currently [human](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9606), [cow](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9913), [pig](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9823), [carp](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=7962)[^carp], [mouse](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=10090), [*E. coli*](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=562)).
- URLs to sequence databases for small and large ribosomal subunits from [SILVA](https://www.arb-silva.de/download/arb-files/).
- Up-to-date links to [VirusHostDB](https://www.genome.jp/virushostdb) and [NCBI Taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/).

[^carp]: The carp genome is included as a last-ditch attempt to [capture any remaining Illumina adapter sequences](https://dgg32.medium.com/carp-in-the-soil-1168818d2191) before moving on to HV identification. I'm not especially confident this is helpful.

Given these inputs, the index workflow:
- Generates a TSV of viral taxa, incorporating taxonomic information from NCBI, and annotates each taxon with infection status for each host taxon of interest (using Virus-Host-DB).
- Makes Bowtie2 indices from (1) all host-infecting viral genomes in Genbank[^genbank], (2) the human genome, (3) common non-human contaminants, plus BBMap indices for (2) and (3).
- Downloads and extracts local copies of (1) the BLAST nt database, (2) the specified Kraken2 DB, (3) the SILVA rRNA reference files.

[^genbank]: Excluding transgenic, contaminated, or erroneous sequences, which are excluded according to a list of sequence ID patterns specified in the config file.

Run the index workflow by setting `mode = "index"` in the relevant config file. For more information, see `workflows/index.nf` and the associated subworkflows and modules.

### Run workflow

#### Preprocessing phase

The run workflow begins by undergoing cleaning by [FASTP](https://github.com/OpenGene/fastp), which both screens for adapters and trims low-quality and low-complexity sequences. The output of FASTP is then passed on to the other parts of the run workflow.

#### Viral identification phase

The goal of this phase is to sensitively, specifically, and efficiently identify reads arising from host-infecting viruses. It uses a multi-step process designed to keep computational costs low while minimizing false-positive and false-negative errors.

1. To begin with, the cleaned reads produced by the preprocessing phase are screened against a database of host-infecting viral genomes generated from Genbank by the index workflow. This initial screen is performed using [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), which flags any read that contains at least three 21-mers matching any host-infecting viral genome. The purpose of this initial screen is to rapidly and sensitively identify putative host-infecting viral reads while discarding the vast majority of non-HV reads, reducing the cost associated with the rest of this phase.
2. Surviving reads undergo additional adapter trimming with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) & [Trimmomatic](https://github.com/timflutre/trimmomatic) to remove any residual adapter contamination that might lead to false positive results.
3. Next, reads are aligned to the previously-mentioned database of HV genomes with [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) using quite permissive parameters designed to capture as many putative HV reads as possible. The SAM and FASTQ files are processed to generate new read files containing any read pair for which at least one read matches the HV database.
4. The output of the previous step is passed to a further filtering step, in which reads matching a series of common contaminant sequences are removed. This is done by aligning surviving reads to these contaminants using both Bowtie2 and [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) in series[^filter]. Contaminants to be screened against include reference genomes from human, cow, pig, carp, mouse and *E. coli*, as well as various genetic engineering vectors.
5. Surviving read pairs are deduplicated with [Clumpify](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/)[^dedup] and merged into single sequences through a combination of [BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)[^merge]. Read pairs that fail to merge with BBMerge are concatenated end-to-end with an intervening "N" base.
6. Deduplicated and merged reads are passed to [Kraken2](https://ccb.jhu.edu/software/kraken2/) for taxonomic assignment, using the reference database obtained in the index workflow. We record whether each read was (1) assigned to a host-infecting virus taxon with Kraken, (2) assigned to a non-HV taxon with Kraken, or (3) not assigned to any taxon. All reads in category (2) are filtered out.
7. Finally, reads are assigned a final HV status if they:
    - Are classified as HV by both Bowtie2 and Kraken2; or
    - Are unassigned by Kraken and align to an HV taxon with Bowtie2 with an alignment score above a user-specifed threshold[^threshold].

Following HV status assignment, information for each read pair is processed into a TSV file available in the output directory as `virus_hits_db.tsv.gz`. Finally, the number of read pairs mapping to each detected HV taxon is counted and output as `virus_clade_counts.tsv.gz` for downstream use.

[^filter]: We've found in past investigations that the two aligners detect different contaminant sequences, and aligning against both is more effective at avoiding false positives than either in isolation.
[^dedup]: Which identifies and collapses groups of reads that are identical modulo some specified error rate.
[^merge]: Which aligns and merges read pairs with significant overlap. 
[^threshold]: Specifically, Kraken-unassigned read pairs are classed as HV if, for either read in the pair, S/ln(L) >= T, where S is the best-match Bowtie2 alignment score for that read, L is the length of the read, and T is the value of `params.bt2_score_threshold` specified in the config file.

#### Taxonomic profiling phase

The goal of this phase is to give an overview of the taxonomic composition of the cleaned reads from the preprocessing phase. In particular, it gives an estimate of (i) the fraction of ribosomal reads in the dataset, (ii) the taxonomic breakdown of the dataset at the domain level[^eukarya], and (iii) more detailed abundance estimates for lower-level taxa.

To maintain efficiency, the reads from the preprocessing phase are first subset, by default to 1 million reads per sample. These subset reads are then separated into ribosomal and non-ribosomal read groups using BBDuk, by searching for ribosomal k-mers from the SILVA database obtained during the index workflow. Both sets of read pairs are then merged into single sequences using BBMerge as described above, then are passed to Kraken2 for taxonomic assignment. The output of Kraken2 is then passed to [Bracken](https://ccb.jhu.edu/software/bracken/) for correction and summarization, and the outputs of both Kraken and Bracken are merged across samples to produce single output files.

[^eukarya]: As human is the only eukaryotic genome included in the Standard reference database for Kraken2, all sequences assigned to that domain can be assigned to *Homo sapiens*.

#### BLAST validation phase

To evaluate the performance of the process described in the viral identification phase, it's useful to get some ground-truth information about whether the HV assignments made in that subworkflow are correct. To do this, we use [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to align the putative HV reads output by the previous phase against the core_nt database, then process the output to check whether each sequence had a high-scoring alignment to at least one HV sequence. For computational tractability, this can be performed on only a subset of surviving HV reads (specified by `params.blast_hv_fraction`)[^blast].

[^blast]: Setting `params.blast_hv_fraction` to 0 skips this step altogether.

Currently, this subworkflow performs the BLAST alignment, filters & summarizes the output to a manageable size, and saves the result to the output directory as `blast_hits_paired.tsv.gz`. Actual assessment of HV status and performance evaluation is currently performed external to this workflow.

#### QC and output phase

In this phase, the raw and cleaned reads from the preprocessing phase are analyzed with [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/). This allows us to assess important metrics (e.g. read quality, read numbers, adapter content) in the raw input, and how these metrics are affected by cleaning. Relevant metrics are then extracted from MultiQC outputs into easy-to-parse TSV files[^tsvs] for downstream processing.

[^tsvs]: Specifically, `qc_basic_stats.tsv.gz`, `qc_adapter_stats.tsv.gz`, `qc_quality_base_stats.tsv.gz` and `qc_quality_sequence_stats.tsv.gz`.

### Pipeline outputs

If the pipeline runs to completion, the following output files are expected.

#### Index workflow

1. `output/input/index-params.json`: JSON file giving all the parameters passed to the pipeline (useful for trying to reproduce someone else's results).
2. `output/logging/pipeline-version.txt`: Version of the pipeline with which index directory was created.
3. `output/logging/time.txt`: Start time of index workflow run.
4. `output/logging/trace.txt`: Nextflow trace file containing logging information for each process performed during the workflow run.
5. `output/results/core_nt`: Directory containing extracted BLAST database files for BLASTing against core_nt.
6. `output/results/bt2-virus-index`: Directory containing Bowtie2 index for host-infecting viral genomes.
7. `output/results/bt2-human-index`: Directory containing Bowtie2 index for the human genome.
8. `output/results/bt2-other-index`: Directory containing Bowtie2 index for other contaminant sequences.
9. `output/results/bbm-human-index`: Directory containing BBMap index for the human genome.
10. `output/results/bbm-other-index`: Directory containing BBMap index for other contaminant sequences.
11. `output/results/kraken_db`: Directory containing Kraken2 reference database (by default, the most recent version of Standard).
12. `output/results/virus-genomes-filtered.fasta.gz`: FASTA file containg host-infecting viral genomes downloaded from viral Genbank (filtered to remove transgenic, contaminated, or erroneous sequences).
13. `virus-genome-metadata-gid.tsv.gz`: Genome metadata file generated during download of HV genomes from viral Genbank, annotated additionally with Genome IDs used by Bowtie2 (allowing mapping between genome ID and taxid).
14. `output/results/ribo-ref-concat.fasta.gz`: Reference database of ribosomal LSU and SSU sequences from SILVA.
15. `output/results/taxonomy-nodes.dmp`: Taxonomy dump file from NCBI mapping between taxids and their parents in the NCBI taxonomy tree structure.
16. `output/results/taxonomy-names.dmp`: Taxonomy dump file from NCBI mapping between taxids and taxon names.
17. `output/results/total-virus-db-annotated.tsv.gz`: Database generated from NCBI taxonomy and Virus-Host-DB giving taxonomy and host-infection information for each viral taxon.

#### Run workflow

1. `output/input`: Directory containing saved input information (useful for trying to reproduce someone else's results)
    1. `adapters.fasta`: FASTA file of adapter sequences used for adapter screening.
    2. `run-params.json`: JSON file giving all the parameters passed to the pipeline.
    3. `index-params.json`: JSON file giving parameters used to generate index directory (`params.ref_dir`).
    4. `samplesheet.csv`: Copy of the samplesheet file used to configure the pipeline (specified by `params.sample_sheet`).
2. `output/logging`: Log files containing meta-level information about the pipeline run itself.
    1. `pipeline-version.txt`: Version of the pipeline used for the run.
    2. `time.txt`: Start time of the run.
    3. `trace.txt`: Tab delimited log of all the information for each task run in the pipeline including runtime, memory usage, exit status, etc. Can be used to create an execution timeline using the the script `bin/plot-timeline-script.R` after the pipeline has finished running. More information regarding the trace file format can be found [here](https://www.nextflow.io/docs/latest/reports.html#trace-file).
    4. `pipeline-version-index.txt`: Version of pipeline used to generate index directory (`params.ref_dir`).
3. `output/intermediates`: Intermediate files produced by key stages in the run workflow, saved for nonstandard downstream analysis.
    1. `reads/cleaned`: Directory containing paired FASTQ files for cleaned reads (i.e. the output of the preprocessing phase described above).
    2. `reads/raw_viral`: Directory containing raw reads corresponding to those reads that survive initial viral screening with BBDuk.
4. `output/results`: Directory containing processed results files for standard downstream analysis.
    1. `qc_basic_stats.tsv.gz`: Summary statistics for each sample at each stage of the preprocessing phase (`stage`), including:
        - GC content (`percent GC`);
        - Average read length (`mean_seq_len`);
        - Number of read pairs (`n_read pairs`);
        - Approximate number of base pairs in reads (`n_bases_approx`);
        - Percent duplicates as measured by FASTQC (`percent_duplicates`);
        - Pass/fail scores for each test conducted by FASTQC.
    2. `qc_adapter_stats.tsv.gz`: Adapter statistics calculated by FASTQC for each sample and preprocessing stage, given as a percentage of reads containing adapter content (`pc_adapters`) at each position along the read (`position`) for each adapter detected (`adapter`) for each read in the read pair (`read_pair`).
    3. `qc_quality_base_stats.tsv.gz`: Per-base read-quality statistics calculated by FASTQC for each sample and preprocessing stage, given as the mean Phred score (`mean_phred_score`) at each position along the read (`position`) for each read in the read pair (`read_pair`).
    4. `qc_quality_sequence_stats.tsv.gz`: Per-sequence read-quality statistics calculated by FASTQC for each sample and preprocessing stage, given as the number of reads (`n_sequences`) with a given mean Phred score (`mean_phred_score`) for each read in the read pair (`read_pair`).
    5. `virus_hits_db.tsv.gz`: TSV output by the viral identification phase, giving information about each read pair assigned to a host-infecting virus.
    6. `virus_clade_counts.tsv.gz`: Summary of the previous file giving the number of HV read pairs mapped to each viral taxon. Includes both read pairs mapped directly to that taxon (`n_reads_direct`) and to that taxon plus all descendent taxa (`n_reads_clade`).
    3. `blast_hits_paired.tsv.gz`: Summarized BLASTN output for putative HV read pairs, giving, for each read pair and subject taxid:
        - The number of reads in the read pair with high-scoring matches to that taxid (`n_reads`).
        - The best bitscores of alignments to that taxid for each matching read (`bitscore_max` and `bitscore_min`)[^bitscore].
    8. `kraken_reports_merged.tsv.gz`: Kraken output reports in TSV format, labeled by sample and ribosomal status.
    9. `bracken_reports_merged.tsv.gz`: Bracken output reports in TSV format, labeled by sample and ribosomal status.

[^bitscore]: If only one read aligns to the target, these two fields will be identical. If not, they will give the higher and lower of the best bitscores for the two reads in the pair..

## Using the workflow

### Profiles and modes

The pipeline can be run in two modes, "index" and "run", each of which calls the corresponding workflow from the `workflows` directory.

For both modes, the pipeline can be run in multiple ways by modifying various configuration variables specified in `configs/profiles.config`. Currently, three profiles are implemented, all of which assume the workflow is being launched from an AWS EC2 instance:

- The `batch` profile is the default and attempts to run the pipeline with AWS Batch. This is the quickest and most efficient way to run the pipeline, but requires significant additional setup not described in this repo. To set up AWS Batch for this pipeline, follow the instructions [here](https://data.securebio.org/wills-public-notebook/notebooks/2024-06-11_batch.html) (steps 1-3), then modify your config file to point `process.queue` to the name of your Batch job queue.
- The `ec2_local` profile attempts to run the whole workflow locally on your EC2 instance, storing intermediate and outflow files on instance-linked block storage. This is simple and can be relatively fast, but is bottlenecked by your instance's CPU and memory allocations; in particular, if you don't use an instance with very high memory, the pipeline is likely to fail when loading the Kraken2 reference DB.
- The `ec2_s3` profile runs the pipeline on your EC2 instance, but attempts to read and write files to a specified S3 directory. This avoids problems caused by insufficient storage on your EC2 instance, but (1) is significantly slower and (2) is still constrained by your instance's memory allocation.

To run the pipeline with a specified profile, run `nextflow run PATH_TO_REPO_DIR -profile PROFILE_NAME -resume`. Calling the pipeline without specifying a profile will run the `batch` profile by default. Future example commands in this README will assume you are using Batch; if you want to instead use a different profile, you'll need to modify the commands accordingly.

> [!TIP]
> It's highly recommended that you always run `nextflow run` with the `-resume` option enabled. It doesn't do any harm if you haven't run a workflow before, and getting into the habit will help you avoid much sadness when you want to resume it without rerunning all your jobs.

### Installation & setup

#### 1. Install dependencies

To run this workflow with full functionality, you need access to the following dependencies:

1. **SDKMAN!:** To install the SDKMAN! Java SDK manager, follow the installation instructions available [here](https://sdkman.io/install).
2. **Nextflow:** To install the workflow management framework, follow the installation instructions available [here](https://www.nextflow.io/docs/latest/getstarted.html), beginning by installing a recommended Java distribution through SDKMAN!. Pipeline version 2.5.0+ requires Nextflow version 24.10.0+.
2. **Docker:** To install Docker Engine for command-line use, follow the installation instructions available [here](https://docs.docker.com/engine/install/) (or [here](https://docs.aws.amazon.com/serverless-application-model/latest/developerguide/install-docker.html) for installation on an AWS EC2 instance).
3. **AWS CLI:** If not already installed, install the AWS CLI by following the instructions available [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
4. **Git:** To install the Git version control tool, follow the installation instructions available [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

#### 2. Configure AWS & Docker

To run the workflow using AWS S3 for the working and output directories, you first need to [configure AWS access](https://www.nextflow.io/docs/latest/aws.html). To do this, you need to create a file at `~/.aws/config` or `~/.aws/credentials` specifying your access key ID and secret access key, e.g.

`~/.aws/config`:
```
[default]
region = us-east-1
output = table
tcp_keepalive = true
aws_access_key_id = <ACCESS_KEY_ID>
aws_secret_access_key = <SECRET_ACCESS_KEY>
```

`~/.aws/credentials`:
```
[default]
aws_access_key_id = <ACCESS_KEY_ID>
aws_secret_access_key = <SECRET_ACCESS_KEY>
```

If you encounter `AccessDenied` errors after doing this, you may also need to export these keys as environment variables before running Nextflow:

```
eval "$(aws configure export-credentials --format env)"
```

Next, you need to make sure your user is configured to use Docker. To do this, create the `docker` user group and add your current user to it:

```
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

#### 3. Clone this repository

Clone this repo into a new directory as normal.

#### 4. Run index/reference workflow

> [!TIP]
> If someone else in your organization already uses this pipeline, it's likely they've already run the index workflow and generated an output directory. If this is the case, you can reduce costs and increase reproducibility by using theirs instead of generating your own. If you want to do this, skip this step, and edit `configs/run.config` such that `params.ref_dir` points to `INDEX_DIR/output`.

Create a new directory outside the repo directory and copy over the index workflow config file as `nextflow.config` in that directory:

```
mkdir index
cd index
cp REPO_DIR/configs/index.config nextflow.config
```

Next, edit `nextflow.config` such that `params.base_dir` points to the directory (likely on S3) where you want to store your index files. (You shouldn't need to modify anything else about the config file at this stage.)

Next, call `nextflow run` pointing at the repo directory (NB: *not* `main.nf` or any other workflow file; pointing to the directory will cause Nextflow to automatically run `main.nf` from that directory.):

```
nextflow run PATH_TO_REPO_DIR -resume
```

Wait for the workflow to run to completion; this is likely to take several hours at least.

### Testing & validation

To confirm that the pipeline works in your hands, we provide a small test dataset (`test/raw`) to run through the run workflow. This can be used to test any of the pipeline profiles described above.

If your EC2 instance has the resources to handle it, the simplest way to start using the pipeline is to run the test data through it locally on that instance (i.e. without using S3). To do this:

1. Navigate to the `test` directory.
2. Edit `nextflow.config` to set `params.ref_dir` to the index directory you chose or created above (specifically `PATH_TO_REF_DIR/output`).
3. Still within the `test` directory, run `nextflow run -profile ec2_local .. -resume`.
4. Wait for the workflow to finish. Inspect the `output` directory to view the processed output files.

If this is successful, the next level of complexity is to run the workflow with a working directory on S3. To do this:

1. Within the `test` directory, edit `nextflow.config` to set `params.base_dir` to the S3 directory of your choice.
2. Still within that directory, run `nextflow run -profile ec2_s3 .. -resume`.
3. Wait for the workflow to finish, and inspect the output on S3.

Finally, you can run the test dataset through the pipeline on AWS Batch. To do this, configure Batch as described [here](https://data.securebio.org/wills-public-notebook/notebooks/2024-06-11_batch.html) (steps 1-3), then:

1. Within the `test` directory, edit `nextflow.config` to set `params.base_dir` to a different S3 directory of your choice and `process.queue` to the name of your Batch job queue.
2. Still within that directory, run `nextflow run -profile batch .. -resume` (or simply `nextflow run .. -resume`).
3. Wait for the workflow to finish, and inspect the output on S3.

### Running on new data

To run the workflow on another dataset, you need:

1. Accessible raw data files in Gzipped FASTQ format, named appropriately.
2. A sample sheet file specifying the samples, along with paths to the forward and reverse read files for each sample.
3. A config file in a clean launch directory, pointing to:
    - The directory containing the raw data (`params.raw_dir`).
    - The base directory in which to put the working and output directories (`params.base_dir`).
    - The directory containing the outputs of the reference workflow (`params.ref_dir`).
    - The sample sheet (`params.sample_sheet`).
    - Various other parameter values.

> [!NOTE]
> The samplesheet must have the following format for each row:
> - First column: Sample ID
> - Second column: Path to FASTQ file 1 which should be the forward read for this sample
> - Third column: Path to FASTQ file 2 which should be the reverse read for this sample
> 
> The easiest way to get this file is by using the `generate_samplesheet.sh` script. As input, this script takes a path to raw FASTQ files (`dir_path`), and forward (`forward_suffix`) and reverse (`reverse_suffix`) read suffixes, both of which support regex, and an optional output path (`output_path`). Those using data from s3 should make sure to pass the `s3` parameter. Those who would like to group samples by some metadata can pass a path to a CSV file containing a header column named `sample,group`, where each row gives the sample name and the group to group by (`group_file`), edit the samplesheet manually after generation (since manually editing the samplesheet will be easier when the groups CSV isn't readily available), or provide the --group_across_illumina_lanes option if a single library was run across a single Illumina flowcell. As output, the script generates a CSV file named (`samplesheet.csv` by default), which can be used as input for the pipeline.


If running on Batch, a good process for starting the pipeline on a new dataset is as follows:

1. Process the raw data to have appropriate filenames (see above) and deposit it in an accessible S3 directory.
2. Create a clean launch directory and copy `configs/run.config` to a file named `nextflow.config` in that directory.
3. Create a library metadata file in that launch directory, specifying library/sample mappings and any other metadata (see above).
4. Edit `nextflow.config` to specify each item in `params` as appropriate, as well as setting `process.queue` to the appropriate Batch queue.
5. Run `nextflow run PATH_TO_REPO_DIR -resume`.
6. Navigate to `{params.base_dir}/output` to view and download output files.

# Troubleshooting

When attempting to run a released version of the pipeline, the most common sources of errors are AWS permission issues. Before debugging a persistent error in-depth, make sure that you have all the permissions specified in Step 0 of [our Batch workflow guide](https://data.securebio.org/wills-public-notebook/notebooks/2024-06-11_batch.html). Next, make sure Nextflow has access to your AWS credentials, such as by running `eval "$(aws configure export-credentials --format env)"`.

Another common issue is for processes to fail with some variation of the following Docker-related error:

```
docker: failed to register layer: write /usr/lib/jvm/java-11-openjdk-amd64/lib/modules: **no space left on device**.
```

This is a fairly well-known problem, which can arise even when there is substantial free storage space accessible to your EC2 instance. Following the steps recommended [here](https://www.baeldung.com/linux/docker-fix-no-space-error) or [here](https://forums.docker.com/t/docker-no-space-left-on-device/69205) typically resolves the issue, either by deleting Docker assets to free up space (e.g. via `docker system prune --all --force`) or by giving Docker more space.

We will add more common failure modes here as they get reported.
