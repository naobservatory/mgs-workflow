# Metagenomic Sequencing Analysis Workflow

This Nextflow workflow is designed to analyze metagenomic sequencing data to characterize overall taxonomic composition and identify and count human-infecting viruses.

## 1. Installation & validation

### 1.1. Install dependencies

To run this workflow with full functionality, you need access to the following dependencies:

1. **SDKMAN!:** To install the SDKMAN! Java SDK manager, follow the installation instructions available [here](https://sdkman.io/install).
2. **Nextflow:** To install the workflow management framework, follow the installation instructions available [here](https://www.nextflow.io/docs/latest/getstarted.html), beginning by installing a recommended Java distribution through SDKMAN!.
2. **Docker:** To install Docker Engine for command-line use, follow the installation instructions available [here](https://docs.docker.com/engine/install/) (or [here](https://docs.aws.amazon.com/serverless-application-model/latest/developerguide/install-docker.html) for installation on an AWS EC2 instance).
3. **AWS CLI:** If not already installed, install the AWS CLI by following the instructions available [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
4. **Git:** To install the Git version control tool, follow the installation instructions available [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

### 1.2. Configure AWS & Docker

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
export AWS_ACCESS_KEY_ID=<ACCESS_KEY_ID>
export AWS_SECRET_ACCESS_KEY=<SECRET_ACCESS_KEY>
```

Next, you need to make sure your user is configured to use Docker. To do this, create the `docker` user group and add your current user to it:

```
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

### 1.3. Run index/reference workflow

Clone this repo into a new directory, then run the index workflow to generate and configure reference and index files. You should only need to run this once for any given S3 bucket. We recommend creating a new directory outside the repo directory and copying/linking the index workflow, config file, and ref directory over to that directory before running the workflow:

```
mkdir index
cd index
ln -s <REPO_DIR>/workflows/index.nf workflow.nf
cp <REPO_DIR>/configs/index.config nextflow.config
cp -r <REPO_DIR>/ref ref
```

Edit the config file to specify the destination S3 bucket, then run the workflow from within the index directory:

```
nextflow workflow.nf -resume
```

Wait for the workflow to run to completion; depending on your computational resources, this could take several hours.

### 1.4. Run local test workflow

To validate functionality and help debug issues, the main workflow comes packaged with a small test dataset and associated configuration file. After running the index/reference workflow, the next step in setting up the workflow for the first time is to run it on this local test dataset.

To run the local test, first edit `<REPO_DIR>/configs/main.config` to point `ref_dir` to the output folder of the index workflow you ran in 1.3. Then navigate to the local test directory and run Nextflow on the main workflow via the pre-existing symbolic link:

```
cd <REPO_DIR>/test-local
nextflow workflow.nf -resume
```

Wait for the workflow to complete successfully.

### 1.5. Run remote test workflow

Finally, confirm complete main workflow functionality by running it on the same test dataset with working and output directories on S3.

```
mkdir test-remote
cd test-remote
ln -s <REPO_DIR>/workflows/main.nf workflow.nf
cp <REPO_DIR>/configs/main.config nextflow.config
cp -r <REPO_DIR>/ref ref
cp -r <REPO_DIR>/scripts scripts
aws s3 cp --recursive <REPO_DIR>/test-local/raw s3://<S3_BUCKET>/<PROJECT_DIR>/raw
```

Next, edit the copied config file to point to your preferred S3 path, as well as updating local paths to the new project directory, by setting the following entries:

```
// Params
s3_dir = "s3://<S3_BUCKET>/<PROJECT_DIR>"
raw_dir = "${params.s3_dir}/raw"
pub_dir = "${params.s3_dir}/output"
script_dir = "${projectDir}/scripts"
library_tab = "${projectDir}/ref/libraries.csv"
adapters = "${projectDir}/ref/adapters.fa"

// Fusion
enabled = true

// Other
workDir = "${params.s3_dir}/work"
```

Then, run Nextflow from the new project directory as normal:

```
nextflow workflow.nf -resume
```

## 2. Running on new data

To run the workflow on new data files, you must first make sure you have access to them from the directory where you intend to run the workflow: for example, via a local path or an accessible S3 bucket. Assuming this is the case, you can run the workflow by following these steps:

### 2.1. Create a new workflow directory

While it isn't strictly necessary, especially when using an S3 working directory, we advise running every instance of the workflow in a new working directory to minimize the possibility of overwriting or otherwise confusing important files.

### 2.2. Prepare library metadata

Next, create a library metadata file specifying the mapping between library files and sample IDs. This should be a CSV file with `library` and `sample` columns:

```
library,sample
D23-14114-1,1A
D23-14114-2,1A
D23-14115-1,2A
D23-14115-2,2A
D23-14116-1,6A
D23-14116-2,6A
...,...
```

You can specify other metadata variables here as well, but it won't currently do anything.

### 2.3. Copy & link workflow files

Copy config and script files from the main repo directory to your run directory, and copy or link the main workflow file:

```
cd <RUN_DIR>
ln -s <REPO_DIR>/workflows/main.nf workflow.nf
cp <REPO_DIR>/configs/main.config nextflow.config
cp -r <REPO_DIR>/ref ref
cp -r <REPO_DIR>/scripts scripts
```

### 2.4. Edit Nextflow config file for test run

When starting to analyze a new dataset, we recommend starting with a test run to (1) confirm workflow functionality and quickly spot any issues, and (2) generate the adapter sequence file for the main run.

To do this, edit `nextflow.config` as follows:

```
// Params
truncate_reads = true
n_reads_trunc = 25000
s3_dir = "s3://<S3_BUCKET>/<PROJECT_DIR>"
raw_dir = "${params.s3_dir}/raw"
pub_dir = "${params.s3_dir}/output"
script_dir = "${projectDir}/scripts"
library_tab = "<PATH_TO_LIBRARY_CSV>"
adapters = "<PATH_TO_ADAPTER_FASTA>"

// Fusion
enabled = true

// Other
workDir = "${params.s3_dir}/work"
```
### 2.5. Execute test run

Execute the test run with the following command:

```
nextflow workflow.nf -resume
```

Wait for completion.

### 2.6. Extract adaptors

If you're confident your adapter sequences are already included in your adapter file from 2.4., you can skip this step.

Otherwise, when the dry run is concluded, copy the inferred adapter file from `<PUB_DIR>/output/results/adapters.fasta` to your run directory. Inspect it to see if it contains any inferred adapter sequences that are missing from your initial adapter file. If it does, add those adapter sequences to your adapter file for the main run.

### 2.7. Edit Nextflow config file for main run

To modify your Nextflow config file for the main run, simply update `params.truncate_reads` to `false` and `params.adapters` to point to your updated adapter file.

### 2.8. Execute main run

As before, execute the main Nextflow run with the same command:

```
nextflow workflow.nf -resume
```

Wait for completion.

### 2.9. Analyze output files

The output files from the Nextflow run are copied to `<PUB_DIR>/output/results/`. At the time of writing, they are as follows:

- `adapters.fasta`: FASTA file listing inferred adapter sequences detected during preprocessing.
- `qc_basic_stats.tsv`, `qc_adapter_stats.tsv`, `qc_quality_base_stats.tsv`, `qc_quality_sequence_stats.tsv`: TSVs containing collated QC information generated by FASTQC and MultiQC.
- `taxonomic_composition.tsv`, `bracken_counts.tsv`: TSVs containing high-level taxonomic composition information for each sample.
- `hv_hits_putative_all.tsv`: TSV listing all putative human-viral reads identified by the Bowtie/Kraken pipeline and associated metadata.
- `hv_hits_putative_filtered.tsv`: TSV listing all putative human-viral reads identified by the Bowtie/Kraken pipeline that pass initial permissive filters on taxonomic assignment & alignment score.

An example of how to use these files to analyze a dataset can be found [here](https://data.securebio.org/wills-public-notebook/notebooks/2023-12-19_project-runway-bmc-rna.html).

## 3. Major outstanding issues

This workflow is intended to ultimately replace [this one](https://github.com/naobservatory/mgs-pipeline) as the primary MGS analysis pipeline used by the NAO team. However, there are some important features that aren't yet implemented:
- Most obviously, this pipeline doesn't contain the terminal steps needed to integrate the results with our dashboards (e.g. [here](https://data.securebio.org/mgs-counts/#mr)).
- It also lacks the ability to re-run across multiple datasets/bioprojects at once, which might be desirable once a large amount of data has been processed with the workflow.

Beyond these points, the current version of the pipeline is known to generate a [significant number of false-positive viral sequences](https://data.securebio.org/wills-public-notebook/notebooks/2024-02-08_crits-christoph-2.html), and requires manual inspection of the results before strong conclusions can be drawn. This issue is under active development.
