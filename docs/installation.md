# Installation & setup

This page describes how to install the pipeline and configure your computing environment to run it.

> [!IMPORTANT]
> This pipeline has been developed to run on Linux environments, primarily on AWS EC2 instances. It has not been tested or debugged on MacOS or Windows machines, and is not guaranteed to work on those environments.

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->
- [Installation \& setup](#installation--setup)
      - [1. Install dependencies](#1-install-dependencies)
      - [2. Configure AWS \& Docker](#2-configure-aws--docker)
      - [3. Clone this repository](#3-clone-this-repository)
      - [4. Run index/reference workflow](#4-run-indexreference-workflow)
<!-- TOC end -->

## 1. Install dependencies 

To run this workflow with full functionality, you need access to the following dependencies:

1. **SDKMAN!:** To install the SDKMAN! Java SDK manager, follow the installation instructions available [here](https://sdkman.io/install).
2. **Nextflow:** To install the workflow management framework, follow the installation instructions available [here](https://www.nextflow.io/docs/latest/getstarted.html), beginning by installing a recommended Java distribution through SDKMAN!. Pipeline version 2.5.0+ requires Nextflow version 24.10.0+.
2. **Docker:** To install Docker Engine for command-line use, follow the installation instructions available [here](https://docs.docker.com/engine/install/) (or [here](https://docs.aws.amazon.com/serverless-application-model/latest/developerguide/install-docker.html) for installation on an AWS EC2 instance).
3. **AWS CLI:** If not already installed, install the AWS CLI by following the instructions available [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
4. **Git:** To install the Git version control tool, follow the installation instructions available [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
5. **nf-test**: To install nf-test, follow the install instructions available [here](https://www.nf-test.com/installation/).

## 2. Configure AWS & Docker

To run the workflow using AWS S3 for the working and output directories, you first need to [configure AWS access](https://www.nextflow.io/docs/latest/aws.html). To do this, you need to create a file at `~/.aws/config` **or** `~/.aws/credentials` specifying your access key ID and secret access key, e.g.

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

> [!TIP]
> If you encounter `AccessDenied` errors after doing this, you may also need to export these keys as environment variables before running Nextflow:
>
> ```
> eval "$(aws configure export-credentials --format env)"
> ```


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

Next, call `nextflow run` pointing at the repo directory:

```
nextflow run PATH_TO_REPO_DIR -resume
```

> [!TIP]
> You don't need to point `nextflow run` at `main.nf` any other workflow file; pointing to the directory will cause Nextflow to automatically run `main.nf` from that directory.

Wait for the workflow to run to completion; this is likely to take several hours at least.
