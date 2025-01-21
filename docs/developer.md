# Developer guide
This section is soley for developers of the pipeline. We thank you greatly for your work!

All of our testing is done using [nf-test](https://www.nf-test.com/), which is a testing framework that allows us to run tests locally and on Github Actions. All tests dataset are stored in the `s3://nao-testing` s3 bucket, and must be made public for nf-test (and anyone) to access. 

As of right now we have the following test datasets ([more details here](#available-test-datasets)):

- `gold-standard-test`: A small test dataset to test the run workflow.

### Contributing to the pipeline

When developing new features, always make sure to make a new branch that starts with your name, then is followed by a one or two word description of the feature you're working on by underscores. For example, `will_streaming` or `harmon_remove_trimmomatic`. This branch should be made from the `dev` branch.

>[!CAUTION]
> Do not make pull requests from the `dev` or `main` branches.

By default, all changes are made to local brnaches, pulled into `dev` then after enough `dev` changes, merged into `main`.

#### Pull request requirements

To have a pull request accepted, you must do the following:

1. **Write new tests** for the changes that you make using `nf-test` if those tests don't already exist. At the very least, these tests should check that the new implementation runs to completion, however tests that verify the output of the test dataset are strongly encouraged.
2. **Pass all existing tests locally** by running `nf-test test tests`. *If you make any changes that affect the output of the pipeline, you must list the changes that have occured in the output of the test dataset during your pull request.*
3. **Update the `CHANGELOG.md` file** with the changes that you are making. More information on how to update the `CHANGELOG.md` file can be found [here](./version_schema.md).
4. **Pass automated tests on Github Actions**. These run automatically when you open a pull request.
5. **Reviewed by a maintainer**.

More information on how to run nf-test locally can be found [here](#running-tests-locally-with-nf-test).

### Everything testing

#### Making a test dataset

To make a new test dataset, copy the test dataset to `s3://nao-testing/[name-of-test-dataset]` where `[name-of-test-dataset]` is the name of the test dataset you want to create (request Harmon, or Will for permission to add to the bucket):

```
aws s3 cp /path/to/my_dataset s3://nao-testing/my_dataset/ --acl public-read
```

> [!NOTE]
> Any time you update a test dataset, you must make it public again.

#### Running tests locally with `nf-test`

By default, we use the [gold standard test](#gold-standard-test-s3nao-testinggold-standard-test) dataset to test the run workflow. To run the tests locally, you need to make sure that you have a big enough ec2-instance. We recommend the `m5.xlarge` with at least `32GB` of EBS storage, as this machine closely reflects the VMs on Github Actions. Once you have an instance, run `nf-test test tests`, which will run all tests in the `tests` directory. If you want to run a specific test, you can either reference the name of the test file, or the tag of the test you want to run:

```
nf-test test tests/main.test.nf # Runs all tests in the main.test.nf file
nf-test test --tag run     # Runs the run workflow
```

The intended results for the run workflow can be found in following directory `test-data/gold-standard-results`. Should the `run_output` test fail, you can diff the resulting files of that test, with the files in this folder to find the differences.

> [!NOTE]
> Periodically delete docker images to free up space on your instance. Running the following command will delete all docker images on your system:
> ```
> docker kill $(docker ps -q) 2>/dev/null || true
> docker rm $(docker ps -a -q) 2>/dev/null || true
> docker rmi $(docker images -q) -f 2>/dev/null || true
> docker system prune -af --volumes
> ```

### Available test datasets

#### Gold standard test (`s3://nao-testing/gold-standard-test/`)
This is the default test dataset that we use for testing the run workflow. It is a small dataset that contains XXX reads from the [Yang 2020](https://www.sciencedirect.com/science/article/abs/pii/S0048969720358514?via%3Dihub) study. The code to generate this data can be found [here](https://github.com/naobservatory/generate-test-dataset/tree/harmon-dev-edits).