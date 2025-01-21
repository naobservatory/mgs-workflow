# Usage

To run any of the workflows you must choose a profile, and then have data to run the pipeline on. To help with this process, we provide a small test dataset to run through the `run` workflow.

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->
- [Usage](#usage)
  - [Profiles and modes](#profiles-and-modes)
    - [Compute resource requirements](#compute-resource-requirements)
  - [Running on new data](#running-on-new-data)
  - [Tutorial: Running the `run` workflow on a test dataset](#tutorial-running-the-run-workflow-on-a-test-dataset)
    - [Setup](#setup)
    - [Profile specific instructions](#profile-specific-instructions)
      - [`ec2_local`](#ec2_local)
      - [`ec2_s3`](#ec2_s3)
      - [`batch`](#batch)
<!-- TOC end -->

## Profiles and modes

The pipeline has three main workflows: `INDEX`, `RUN`, and `RUN_VALIDATION`. Each of these calls their corresponding subworkflows, which are located in the `subworkflows` directory.

The pipeline can be run in multiple ways by modifying various configuration variables specified in `configs/profiles.config`. Currently, three profiles are implemented, all of which assume the workflow is being launched from an AWS EC2 instance:

- `batch (default)`:  Most efficient way to run the pipeline
  - This profile is the default and attempts to run the pipeline with AWS Batch. This is the quickest and most efficient way to run the pipeline, but requires significant additional setup not described in this repo. To set up AWS Batch for this pipeline, follow the instructions [here](https://data.securebio.org/wills-public-notebook/notebooks/2024-06-11_batch.html) (steps 1-3), then modify your config file to point `process.queue` to the name of your Batch job queue.
- `ec2_local`: Simple and can be relatively fast, but is bottlenecked by your instance's CPU and memory allocations.
  - This profile attempts to run the whole workflow locally on your EC2 instance, storing intermediate and outflow files on instance-linked block storage. This is simple and can be relatively fast, but is bottlenecked by your instance's CPU and memory allocations; in particular, if you don't use an instance with very high memory, the pipeline is likely to fail when loading the Kraken2 reference DB.
- `ec2_s3`: Avoids storage issues on your EC2 instance, but is still constrained by your instance's memory allocation.
  - This profile runs the pipeline on your EC2 instance, but attempts to read and write files to a specified S3 directory. This avoids problems caused by insufficient storage on your EC2 instance, but (1) is significantly slower and (2) is still constrained by your instance's memory allocation.

To run the pipeline with a specified profile, run `nextflow run PATH_TO_REPO_DIR -profile PROFILE_NAME -resume`. Calling the pipeline without specifying a profile will run the `batch` profile by default. Future example commands in this README will assume you are using Batch; if you want to instead use a different profile, you'll need to modify the commands accordingly.

As of the current pipeline, you must have at least 128GB of memory and 32 cores to run the pipeline if you're running ec2_local or ec2_s3. This is because we use the whole KrakenDB whihc is large (128GB) and some for processes consume 64 cores. Simiarly, if one would like to run BLAST, they must have at least 256GB of memory. 

> [!TIP]
> It's highly recommended that you always run `nextflow run` with the `-resume` option enabled. It doesn't do any harm if you haven't run a workflow before, and getting into the habit will help you avoid much sadness when you want to resume it without rerunning all your jobs.

### Compute resource requirements

To run the pipeline as is you need at least 128GB of memory and 64 cores. This is because we use the whole KrakenDB whihc is large (128GB) and some for processes consume 64 cores. Simiarly, if one would like to run BLAST, they must have at least 256GB of memory. 

To change the compute resources for a process, you can modify the `resources.config` file. This file specifies the compute resources for each process based on the label of the process. For example, to change the compute resources for the `kraken` process, you can add the following to the `resources.config` file:

In the case that you change the resources, you'll need to also change the index.


## Running on new data

To run the workflow on new data, you need:

1. Accessible raw data files in Gzipped FASTQ format, named appropriately.
2. A sample sheet file specifying the samples, along with paths to the forward and reverse read files for each sample. `generate_samplesheet.sh` (see below) can make this for you.
3. A config file in a clean launch directory, pointing to:
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
>
> For example:
> ```
> ../bin/generate_samplesheet.sh \
>    --s3
>    --dir_path s3://nao-restricted/MJ-2024-10-21/raw/ \
>    --forward_suffix _1 \
>    --reverse_suffix _2
> ```

If running on Batch, a good process for starting the pipeline on a new dataset is as follows:

1. Process the raw data to have appropriate filenames (see above) and deposit it in an accessible S3 directory.
2. Create a clean launch directory and copy `configs/run.config` to a file named `nextflow.config` in that directory.
3. Create a sample sheet in that launch directory (see above)
4. Edit `nextflow.config` to specify each item in `params` as appropriate, as well as setting `process.queue` to the appropriate Batch queue.
5. Run `nextflow run PATH_TO_REPO_DIR -resume`.
6. Navigate to `{params.base_dir}/output` to view and download output files.


## Tutorial: Running the `run` workflow on a test dataset

To confirm that the pipeline works in your hands, we provide a small test dataset (`s3://nao-testing/gold-standard-test/raw/`) to run through the `run` workflow. Feel free to use any profile, but we recommend using the `ec2_local` profile as long as [you have the resources](usage.md#compute-resource-requirements) to handle it. 

When running with any profile, there is a shared setup that you need to do, before running the workflow. The shared setup is the same for all profiles, and is described below, followed by profile specific instructions.

### Setup

1. Create a new directory outside the repo directory and copy over the run workflow config file as `nextflow.config` in that directory:

```
mkdir launch
cd launch
cp REPO_DIR/configs/run.config nextflow.config
```

2. Edit `nextflow.config` to set `params.ref_dir` to the index directory you chose or created above (specifically `PATH_TO_REF_DIR/output`)
3. Set the samplesheet path to the test dataset samplesheet `${projectDir}/test-data/samplesheet.csv`

### Profile specific instructions

#### `ec2_local`

4. Within this directory, run `nextflow run -profile ec2_local .. -resume`. Wait for the workflow to finish. 
5. Inspect the `output` directory to view the processed output files.

#### `ec2_s3`

4. Edit `nextflow.config` to set `params.base_dir` to the S3 directory of your choice. 
5. Still within that directory, run `nextflow run -profile ec2_s3 .. -resume`. 
6. Wait for the workflow to finish, and inspect the output on S3.

#### `batch`

4. Edit `nextflow.config` to set `params.base_dir` to the S3 directory of your choice and `process.queue` to the name of your Batch job queue. 
5. Still within that directory, run `nextflow run -profile batch .. -resume` (or simply `nextflow run .. -resume`). 
6. Wait for the workflow to finish, and inspect the output on S3.