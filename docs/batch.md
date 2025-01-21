# AWS BATCH

This page is a guide for running the pipeline on [AWS Batch](https://aws.amazon.com/batch/), a tool which allows you to run Nextflow workflows in a highly parallelized and automated manner. The original source of this guide is [this notebook](https://data.securebio.org/wills-public-notebook/notebooks/2024-06-11_batch.html), however this verison will be updated to reflect the current state of the pipeline. 

## 0. Check your permissions

The most common ways you can run into trouble following these instructions all arise from insufficient AWS permissions. For a smooth experience, make sure you have the following:

1. **AWS Batch Permissions:** You need to have appropriate permissions to view, modify and run Batch compute environments and job queues. The simplest way to do this is to have your administrator add you to the `AWSBatchFullAccess` IAM policy (at SecureBio, this is enabled for all members of the `TeamMembers` group).
2. **S3 Permissions:** You need to have appropriate permissions to read, write and view the S3 bucket where your workflow data is stored. A useful start here is to have your administrator add you to the `AmazonS3FullAccess` IAM policy. After that, you should make sure you have the appropriate permissions to access the specific bucket you're using for your workflow, including the ability to list, read, write, and delete objects[^s3].
3. **EC2 Permissions:** You need to have appropriate permissions to create and modify EC2 launch templates and compute environments. The simplest way to do this is to have your administrator add you to the `AmazonEC2FullAccess` IAM policy.
4. **Instance Role Permissions:** You will need to assign an Instance Role when setting up your Batch compute environment. This role should have at least the following permissions: `AmazonEC2ContainerServiceRole`, `AmazonEC2ContainerServiceforEC2Role`, and `AmazonS3FullAccess`. Make sure you can either set up such a role yourself, or have your administrator do so and point you to the role name. (On SecureBio, this role is called `ecsInstanceRole`.)

[^s3]: In more depth, you need the following actions to be enabled for the bucket in question for your IAM user or role: `s3:ListBucket`, `s3:GetBucketLocation`, `s3:GetObject`, `s3:GetObjectAcl`, `s3:PutObject`, `s3:PutObjectAcl`, `s3:PutObjectTagging`, and `s3:DeleteObject`. If you're using a bucket specific to your user, all this is easier if you first have your administrator enable `s3:GetBucketPolicy` and `s3:PutBucketPolicy` for your user.

## 1. Create an EC2 launch template

First, you need to create an EC2 launch template that specifies the configuration for EC2 instances to be set up through Batch. To do this on the AWS console, Navigate to EC2 on the Services menu, then:

1.  Select "Launch templates" in the left-hand menu.
2.  Click the orange "Create launch template" button.
3.  Enter a launch template name (e.g. "will-batch-template") and optionally a description.
4.  Check the box under "Auto scaling guidance"
5.  Under "Application and OS Images (Amazon Machine Image)", click "Browse more AMIs", then search for "Amazon ECS-Optimized Amazon Linux 2023 (AL2023) x86_64 AMI".
    1.  Under "AWS Marketplace AMIs", select "Amazon ECS-Optimized Amazon Linux 2023 (AL2023) x86_64 AMI" by Amazon Web Services.
    2.  In the popup, select "Subscribe now"
6.  Select an instance type (this isn't hugely important as Batch will modify the instance types provisioned based on the needs of the workflow; I generally select "m5.8xlarge").
7.  Under "Key pair (login)", select "Don't include in launch template"
8.  Under "Network settings", select "Create security group" and follow the default settings.
9.  Now we come to storage. Configuring this correctly is important to avoid IOPS errors!
    1.  The key thing to realize is that, since Batch is spinning up and terminating instances as needed, the usual costs of creating large EBS volumes don't really apply. As such, you can be relatively greedy in provisioning storage for these instances, to minimize the risk of IOPS-related problems with your workflow.
    2.  My recommendation is as follows: under "Storage (volumes)", expand the default "Volume 1 (AMI Root)" volume, then enter the following configuration values[^1]:
        1.  Size: 1000 GiB
        2.  Volume type: gp3
        3.  IOPS: 16000 (maximum for gp3)
        4.  Throughput: 1000
10. Finally, add tags so the cost of your provisioned resources can be tracked more effectively. Add one "Instances" tag (e.g. "will-batch-template-instance") and one "Volumes" tag (e.g. "will-batch-template-volumes").
11. Leave the other settings as default (mostly "Don't include in launch template") and select "Create launch template".

[^1]: If you want even more IOPS, you can provision an io2 volume instead of gp3. However, that's beyond the scope of this guide.

If you want to modify your template after creating it, you can do so by navigating to it in the panel and selecting "Actions" \> "Modify template (create new version)". Be careful to pay attention to which version of the template any dependent resources (e.g. compute environments, see below) are using.

## 2. Set up a Batch compute environment

Next, you need to create a compute environment through which jobs can be allocated to instances. To do this on the AWS console, navigate to Batch on the Services menu, then:

1.  Select "Compute environments" in the left-hand menu
2.  Click the orange "Create" button
3.  Under "Compute environment configuration", select "Amazon Elastic Compute Cloud (Amazon EC2)"[^2]. Then:
    1.  Under "Orchestration Type", select "Managed".
    2.  Enter an environment name (I generally go with something unimaginative like "will-batch-3".
    3.  Set up roles:
        1.  Under "Service role" select "AWSServiceRoleForBatch".
        2.  Under "Instance role" select "ecsInstanceRole", or another role with appropriate permissions (see step 0 above).
        3.  If these roles don't exist or aren't available in the drop-down menu, contact your administrator about setting them up for you.
4.  Under "Tags", navigate to "EC2 tags" and click the "Add tag" button. Then select a tag name that can be used to uniquely identify your use of the workflow (e.g. "mgs-workflow-will"). This is important to let your administrator keep track of how much money you and the wider team are spending on Batch (whose resource consumption is otherwise somewhat inscrutable).
5.  Click the orange "Next" button. Then, under "Instance configuration":
    1.  Under "Use EC2 Spot instances", make sure the "Enable using Spot instances" selector is enabled.
    2.  Under "Maximum % on-demand price", enter a number between 0 and 100. 100 is a good default. Lower numbers will lower costs but increase the chance of unscheduled instance termination, which will require your workflow to re-run jobs.
    3.  Enter integer values under "Minimum vCPUs", "Desired vCPUs" and "Maximum vCPUs". I typically use 0, 0, and 1024.
    4.  Under "Allowed instance types", select "optimal" plus whatever other instance families you want to provision. I typically use optimal, m5, and c6a.
    5.  Under "Additional configuration":
        1.  Specify your EC2 key pair to allow direct ssh'ing into Batch instances (you should very rarely need to do this so you can skip this if you like)
        2.  Under "Launch template" select the launch template you configured previously.
        3.  Under "Launch template version", enter "\$Default" or "\$Latest" (your preference).
6.  Click the orange "Next button", then do so again (i.e. accept defaults for "Network configuration").
    1.  You can configure your own network setup if you like, but that's beyond the scope of this guide.
7.  Review your selections, then click the orange "Create compute environment" button.

[^2]: In the future, I'll investigate running Batch with Fargate for Nextflow workflows. For now, using EC2 gives us greater control over configuration than Fargate, at the cost of additional setup complexity and occasional startup delays.

## 3. Set up a Batch job queue

The last step you need to complete on AWS itself is to set up a job queue that Nextflow can use to submit jobs to Batch. To do this on the AWS console, navigate to Batch on the Services menu, then:

1.  Select "job queues" in the left-hand menu.
2.  Click the orange "Create" button.
3.  Under "Orchestration Type", again select "Amazon Elastic Compute Cloud (Amazon EC2)".
4.  Under "Job queue configuration", enter:
    1.  A queue name (I generally go with something uncreative like "will-batch-queue-nf")
    2.  A job priority (unimportant if you're only using one queue and you're the only one using that queue)
    3.  A connected compute environment (select the environment you set up previously from the dropdown menu)
5.  Click the orange "Create job queue" button.
6.  Success!

## 4. Run Nextflow with Batch

Finally, you need to use all the infrastructure you've just set up to actually run a Nextflow workflow! We recommend using our test dataset to get started. [Click here to see how to run the pipeline on the test dataset](./testing.md#user-guide).