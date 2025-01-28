# Troubleshooting

## Permission issues

When attempting to run a released version of the pipeline, the most common sources of errors are AWS permission issues. Before debugging a persistent error in-depth, make sure you have all the permissions required. When running the pipeline on AWS Batch, the necessary permissions are specified in [our Batch tutorial](./batch.md#step-0-set-up-your-aws-credentials).

Once you've verified you have the required permissions, make sure Nextflow has access to your AWS credentials (if not you may get `AccessDenied` errors):

```
eval "$(aws configure export-credentials --format env)"
```

## Docker image failures

Another common issue is for processes to fail with some variation of the following Docker-related error:

```
docker: failed to register layer: write /usr/lib/jvm/java-11-openjdk-amd64/lib/modules: **no space left on device**.
```

This is a fairly well-known problem, which can arise even when there is substantial free storage space accessible to your EC2 instance. Following the steps recommended [here](https://www.baeldung.com/linux/docker-fix-no-space-error) or [here](https://forums.docker.com/t/docker-no-space-left-on-device/69205) typically resolves the issue, either by deleting Docker assets to free up space (e.g. via `docker system prune --all --force`) or by giving Docker more space.

## Resource constraint errors

Jobs may sometimes fail due to insufficient memory or CPU availability, especially on very large datasets or small instances. To fix this, you can:
- **Increase resource allocations in `configs/resources.config`.** This will alter the resources available to all processes with a given tag (e.g. "small").
- **Increase resource allocation to a specific process.** You can do this by editing the process in the relevant Nextflow file, most likely found at `modules/local/MODULE_NAME/main.nf`.
Note that in some cases it may not be possible to allocate enough resources to meet the needs of a given process, especially on a resource-constrained machine. In this case, you will need to use a smaller reference file (e.g. a smaller Kraken reference DB) or obtain a larger machine.
