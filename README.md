# Nucleic Acid Observatory Viral Metagenomics Pipeline

This Nextflow pipeline is designed to process metagenomic sequencing data, characterize overall taxonomic composition, and identify and quantify reads mapping to human-infecting viruses.

## Pipeline description

### Overview

The pipeline currently consists of two workflows, an **index** workflow and a **run** workflow. The former creates index and reference files that are used by the latter, and is intended to be run first, after which many instantiations of the latter workflow can use the same index output files.

The run workflow then consists of five phases:

1. A **preprocessing phase**, in which input files undergo adapter & quality trimming, deduplication, and ribodepletion.
2. A **taxonomic profiling phase**, in which Kraken2 is used to assess the overall taxonomic composition of the input data and assess how it changes across different steps in the preprocessing phase.
3. A **viral identification phase**, in which a custom multi-step pipeline based around Bowtie2 and Kraken2 is used to sensitively and specifically identify human-infecting virus (HV) reads in the input data for downstream analysis.
4. An optional **BLAST validation phase**, in which putative HV reads from phase 3 are checked against nt to evaluate the sensitivity and specificity of the HV identification process.
5. A final **QC and output phase**, in which FASTQC, MultiQC and other tools are used to assess the quality of the data produced by the pipeline at various steps and produce summarized output data for downstream analysis and visualization.

A slight simplification of the overall process is given by this mildly baffling diagram:

![A flowchart summarizing the pipeline's run workflow.](/readme-workflow-diagram.png)

### Index workflow

The index workflow is run by setting `mode = "index` in the relevant config file, and is intended to be run once (per reasonable length of time) to generate static index and reference files which can then be used by many instantiations of the run workflow. Many of these index/reference files are derived from publicly available reference genomes or other resources, and should thus be updated and re-run periodically as new versions of these become available; however, to keep results comparable across datasets analyzed with the run workflow, this should be done relatively rarely.

Key inputs to the index workflow include:
- A URL linking to a suitable Kraken2 database for taxonomic profiling (typically the [latest release](https://benlangmead.github.io/aws-indexes/k2) of the `k2_standard` database).
- URLS for up-to-date releases of reference genomes for various common contaminant species that can confound the identification of HV reads (currently [human](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9606), [cow](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9913), [pig](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9823), [carp](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=7962), [mouse](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=10090)[^carp], [*E. coli*](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=562)).
- URLs to sequence databases for small and large ribosomal subunits from [SILVA](https://www.arb-silva.de/download/arb-files/).
- Up-to-date links to [VirusHostDB](https://www.genome.jp/virushostdb) and [NCBI Taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/).

[^carp]: The carp genome is included as a last-ditch attempt to [capture any remaining Illumina adapter sequences](https://dgg32.medium.com/carp-in-the-soil-1168818d2191) before moving on to HV identification. I'm not especially confident this is helpful.

Given these inputs, the index workflow:
- Obtains an up-to-date list of human-infecting viruses and their corresponding taxids from VirusHostDB, incorporates information about their taxonomic structure (in particular, their parent taxid) from NCBI, and generates a tap-separated database of human-infecting viral taxa.
- Generates a similar TSV for all viral taxa, including those that don't infect humans.
- Makes Bowtie2 indices from (1) all human-infecting viral genomes in Genbank[^genbank], (2) the human genome, (3) common non-human contaminants, plus BBMap indices for (2) and (3).
- Downloads and extracts local copies of (1) the BLAST nt database, (2) the specified Kraken2 DB, (3) the SILVA rRNA reference files.

[^genbank]: Excluding transgenic, contaminated, or erroneous sequences, which are excluded according to a list of sequence ID patterns specified in the config file.

For more information, see `workflows/index.nf` and the associated subworkflows and modules.

### Run workflow
#### Preprocessing phase

The run workflow begins by concatenating all libraries assigned to the same sample ID together[^concat], then optionally subsampling the concatenated read files down to a specified number of reads[^subsample]. The reads then undergo cleaning by Cutadapt, Trimmomatic, and FASTP, in that order; the first two of these only screen for adapter contamination, while the last both screens for adapters and trims low-quality and low-complexity sequences. Running three different adapter-trimming tools in succession is unusual, but necessary to minimize spurious identification of "human-viral" reads based on adapter contamination of the reference genomes.

[^concat]: This is controlled by the library file specied in the config file; any libraries with the same entry in the `sample` column are concatenated. This is primarily useful in cases where the same library is sequenced multiple times, e.g. to reach some total target depth.
[^subsample]: This is mainly useful for an initial run of the pipeline to confirm everything functions well for a given dataset and sample sheet; it's not recommended for analyses that are actually expected to yield meaningful results.

Cleaned reads then undergo deduplication with Clumpify, a tool from the BBTools suite that identifies and collapses read pairs that are identical modulo some specified error rate. (Importantly, this deduplication phase does not remote reverse-complement duplicates.[^rc]) The deduplicated reads then undergo ribodepletion with BBDuk, which searches for ribosomal k-mers based on a SILVA sequencing database and removes any reads that exhibit rRNA matches above some threshold. This is performed twice in series, once with relatively strict parameters for what sequences qualify as "ribosomal" (which thus removes relatively fewer reads) and again with more inclusive parameters (which thus result in more reads being removed). The output of the initial, more conservative ribodepletion step is passed to the viral identification phase, while the output of the latter, more aggressive step is passed to the taxonomic profiling phase[^ribo].

[^rc]: Unfortunately, I'm not aware of any deduplication tool that (1) works well on paired-read data, (2) can detect and remove duplicates in the presence of sequencing errors, and (3) can handle reverse-complement duplicates. If you know of one, please do let me know!
[^ribo]: The splitting of ribodepletion into these two phases arises from differences in the needs of the two downstream processes they feed into. For viral identification, we want to make sure we aren't losing any HV reads through spurious identification as ribosomal, so we use more conservative parameters. For taxonomic profiling, we're much less concerned about the status of individual reads and more concerned about accurately measuring ribosomal content, so a more aggressive approach is appropriate.

#### Taxonomic profiling phase

In the taxonomic profiling phase, paired-end reads output from the preprocessing phase are merged together into single sequences through a combination of BBMerge (which aligns and merges read pairs with significant overlap) and end-to-end concatenation (with an intervening "N" base, for those read pairs BBMerge is unable to merge). The resulting single sequences are passed to Kraken2 for taxonomic assignment, using the reference database obtained in the index workflow. The output of Kraken2 is then passed to Bracken for correction and summarization, and the outputs of both Kraken and Bracken are merged across samples to produce single output files.

This phase is actually run three times across the pipeline: once on the cleaned reads, once on deduplicated reads, and once on the output of the second ribodepletion step. The latter of these is the primary source of data for overall assessment of the composition of the dataset; the former two are primarily used to assess how deduplication alters the measured taxonomic composition.

#### Viral identification phase
#### BLAST validation phase
#### QC and output phase

TODO

### Pipeline outputs
#### Index workflow
#### Run workflow

TODO

## Using the workflow

TODO: Discuss EC2, S3, and profiles
TODO: Link to Batch guidance

### Installation & setup

#### 1. Install dependencies

To run this workflow with full functionality, you need access to the following dependencies:

1. **SDKMAN!:** To install the SDKMAN! Java SDK manager, follow the installation instructions available [here](https://sdkman.io/install).
2. **Nextflow:** To install the workflow management framework, follow the installation instructions available [here](https://www.nextflow.io/docs/latest/getstarted.html), beginning by installing a recommended Java distribution through SDKMAN!.
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

If you encounter AccessDenied errors after doing this, you may also need to export these keys as environment variables before running Nextflow:

```
eval "$(aws configure export-credentials --format env)
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

#### 4. Run index/reference workflow (if necessary)

> [!TIP]
> If someone else in your organization already uses this pipeline, it's likely they've already run the index workflow and generated an output directory. If this is the case, you can reduce costs and increase reproducibility by using theirs instead of generating your own. If you want to do this, skip this step, and edit `configs/run.config` such that `params.ref_dir` points to `INDEX_DIR/output/results`.

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

> [!TIP]
> It's highly recommended that you always run `nextflow run` with the `-resume` option enabled. It doesn't do any harm if you haven't run a workflow before, and getting into the habit will help you avoid much sadness when you want to resume it without rerunning all your jobs.

> [!NOTE]
> As always, run the above command with `-profile ec2_local` or `-profile ec2_s3` if you don't want to use AWS Batch.

Wait for the workflow to run to completion; this is likely to take several hours at least.

### Testing & validation

To confirm that the pipeline works in your hands, we provide a small test dataset (`test/raw`) to run through the run workflow. This can be used to test any of the pipeline profiles described above.

If your EC2 instance has the resources to handle it, the simplest way to start using the pipeline is to run the test data through it locally on that instance (i.e. without using S3). To do this:

1. Navigate to the `test` directory.
2. Edit `nextflow.config` to set `params.ref_dir` to the index directory you chose or created above (specifically `PATH_TO_REF_DIR/output/results`).
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
2. A library file specifying the relationship between the names of these files (specifically, substrings uniquely identifying each pair of FASTQ files) and the original samples.
3. A sample sheet giving metadata for each sample (can be identical with the library file).
4. A config file in a clean launch directory, pointing to:
    - The directory containing the raw data (`params.raw_dir`).
    - The base directory in which to put the working and output directories (`params.base_dir`).
    - The directory containing the outputs of the reference workflow (`params.ref_dir`).
    - The library file (`params.library_tab`) and sample sheet (`params.sample_tab`).
    - Various other parameter values.

> [!NOTE]
> Currently, the pipeline requires the following of raw data files:
>   - They must be contained in a single directory
>   - Each pair of reads files must have names ending in `_1.fastq.gz` and `_2.fastq.gz`, respectively
>   - Each pair of read files must be uniquely identifiable by a filename substring (specified in the `library` column of `params.library_tab`)

> [!NOTE]
> The library file should specify the mapping between library read files and sample IDs. It must be a CSV file with `library` and `sample` columns, as well as any other metadata columns you feel is appropriate. See `test/libraries.csv` for a minimal example.

If running on Batch, a good process for starting the pipeline on a new dataset is as follows:

1. Process the raw data to have appropriate filenames (see above) and deposit it in an accessible S3 directory.
2. Create a clean launch directory and copy `configs/run.config` to a file named `nextflow.config` in that directory.
3. Create a library metadata file in that launch directory, specifying library/sample mappings and any other metadata (see above).
4. Edit `nextflow.config` to specify each item in `params` as appropriate, as well as setting `process.queue` to the appropriate Batch queue.
5. Run `nextflow run PATH_TO_REPO_DIR -resume`.
6. Navigate to `{params.base_dir}/output` to view and download output files.
