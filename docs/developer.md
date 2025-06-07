# Developer guide
This section is solely for developers of the pipeline. We thank you greatly for your work! It includes guidance on:
- Coding style
- Testing
- GitHub issues
- Pull requests (PRs)
- New releases

The pipeline is written primarily in [Nextflow](https://www.nextflow.io/docs/latest/index.html), with select scripts in Python, Rust, and R. 

## Coding style guide

These guidelines represent best practices to implement in new code, though some legacy code may not yet conform to all conventions.

### Nextflow
- Code organization
    - We use a workflow of workflows organization (`main.nf` -> workflow -> subworkflow -> process). 
    - Aim for one process per module, in a `main.nf` file.
    - Avoid creating duplicate processes. If you need a slight variation on existing behavior, parameterize or otherwise tweak an existing process.
- Documentation
    - Extensive comments are encouraged. 
    - Each workflow or process should begin with a descriptive comment explaining what it does.
    - Each workflow should have a `<workflow_name>.md` document in `docs/`.
- Process conventions (see `modules/local/vsearch/main.nf` for an example of a well-written process that follows these conventions):
    - All processes should have a label specifying needed resources (e.g. `label "small"`). Resources are then specified in `configs/resources.config`.
    - All processes should have a label specifying the Docker container to use (e.g. `label "BBTools"`). Containers are then specified in `configs/containers.config`.
    - Any processes that are used only for testing should have `label "testing"`.
    - All processes should emit their input (for testing validation); use `ln -s` to link the input to the output.
    - Most processes have two output channels, `input` and `output`. If a process emits multiple types of output, use meaningful emit names describing the output types (e.g. `match`, `nomatch`, and `log` from `process BBDUK`).
    - Most, but not all, processes are *labeled* (input is a tuple of a sample name and file path). If input is labeled, output should also be labeled.
- Containers: We preferentially use [Seqera containers](https://seqera.io/containers/), with [Docker Hub](https://hub.docker.com/) as a second choice.
- Naming:
    - Use `lower_snake_case` for variable and channel names.
    - Use `UPPER_SNAKE_CASE` for process and workflow names.
    - Use `camelCase` for module and subworkflow directory names.
    - Use lowercase with hyphens for other file names (config and test data files).
- Performance conventions:
    - Design processes for streaming: avoid loading significant data into memory.
    - Use compressed intermediate files to save disk space.

### Other languages (Python, Rust, R)
- Add non-Nextflow scripts only when necessary; when possible, use existing bioinformatics tools and shell commands rather than creating custom scripts.
- New scripts should be in Rust (preferred) or Python (acceptable if performance is not critical).  We have a few legacy scripts in R but discourage adding new R scripts.
- Organization:
    - Python and R scripts for a module go in `resources/usr/bin/`
    - Rust source code for a module goes in `src` and the binary (compiled with `cargo build --release`) goes in `resources/usr/bin`
        - Note that additional files created by Cargo (`Cargo.lock`, `Cargo.toml`) are not stored in the repo.
- Rely on the standard library as much as possible. (E.g. for Python, avoid use of third-party libraries like `pandas`.)
- Include proper error handling and logging
- Performance conventions:
    - Always process large files line-by-line.
    - Support compressed file formats (.gz, .bz2). 
- Python style: 
    - Loosely follow PEP 8 conventions.
    - Type hints are encouraged but not currently required.
    - Linting is encouraged (our go-to tool is `ruff`), but not currently required. 
    
## Testing

All tests use the [nf-test](https://www.nf-test.com/) framework. (We do not currently have any unit tests of Python/Rust/R scripts.)

### Organization of tests 

Tests themselves are organized in the `tests/` directory following the same structure as the main code. 
```
tests/
├── modules/local/
│   └── toolName/
│       └── main.nf.test
├── subworkflows/local/
│   └── workflowName/
│       └── main.nf.test
└── workflows/
    └── workflow_name.nf.test
```

We also have `.github/workflows/end-to-end.yml` which specifies which tests run automatically on PRs through Github Actions. 

Config files for tests are organized as follows:
- `tests/nextflow.config` is the default config file for `nf-test` and specifies resources for each process based on labels.
- `tests/config/` has config files used by each workflow test.
- `configs/index-for-run-test.config` (NOT in the `tests/` directory) is used to run the `INDEX` workflow (in non-test mode) to generate the index used for `RUN` workflow tests.

Small (uncompressed) test data files are in `test-data/`; larger test datasets are in S3:
- We have some small "toy" test data files in `test-data/toy-data`.
- Larger public test datasets are stored in `s3://nao-testing`. 
- Results of workflow runs on the test datasets from S3 are in the repo in `test-data/<dataset>-results-<workflow>`. 
   
Test datasets are described in more detail in the "Test datasets" section below.

### Writing tests

[Documentation for nf-test](https://www.nf-test.com/docs/getting-started/) and Nextflow's [nf-test training](https://training.nextflow.io/2.1.3/side_quests/nf-test/) are helpful resources.

Here are some guidelines:
- At a minimum, write a test that the process completes successfully on valid input data (`assert process.success`).
- When relevant, processes should be tested on both single- and paired-end data.
- Use descriptive test names and appropriate tags.
- When practical:
     - Validate that output files exist and have expected properties.
     - Include both positive tests (expected success) and edge case/error tests.

`tests/modules/local/fastqc/main.nf.test` is an example of bare-minimum tests for a process, and `tests/modules/local/vsearch/main.nf.test` is an example of really good, comprehensive testing!

### Test datasets 

Large test datasets are in s3://nao-testing (publicly available). 

#### Gold standard test (`s3://nao-testing/gold-standard-test/`)
This is the default test dataset that we use for testing the `RUN` and `DOWNSTREAM` workflows on short-read data. It is a small dataset that contains 165 reads from the [Yang 2020](https://www.sciencedirect.com/science/article/abs/pii/S0048969720358514?via%3Dihub) study. The code to generate this data can be found [here](https://github.com/naobservatory/generate-test-dataset/blob/main/).

#### ONT wastewater test dataset (`s3://nao-testing/ont-ww-test/`)
This is the default, small test dataset that we use for testing the `RUN` workflow on ONT (long-read) data.

### Adding new test data

Small "toy" data files (uncompressed, generally ~1KB or less) may be added freely to the repo in `test-data/toy-data`.

To make a new test dataset on S3, copy the test dataset to `s3://nao-testing/<name-of-test-dataset>`. A pipeline maintainer (e.g. willbradshaw or harmonbhasin) can give you permission to add to the bucket.

```
aws s3 cp /path/to/my_dataset s3://nao-testing/my_dataset/ --acl public-read
```

> [!NOTE]
> Any time you update a test dataset, you must make it public again.

### Running tests

To run the tests locally, you need to make sure that you have a powerful enough compute instance (at least 4 cores, 14GB of RAM, and 32GB of storage). On AWS EC2, we recommend the `m5.xlarge`. Note that you may want a more powerful instance when running tests in parallel (as described below).

To run specific tests, you can specify the tests by filename or by tag. Individual tests generally complete quickly (seconds to minutes):
```
nf-test test tests/main.test.nf # Runs all tests in the main.test.nf file
nf-test test --tag run # Runs test(s) with the "run" tag; this is the end-to-end test of the RUN workflow on short-read data
```

To run all tests in the `tests` directory, use the command `nf-test test tests`. 
- Running this full suite of local tests takes hours; we recommend using the script `bin/run_parallel_test` to parallelize. 
- Running the full test suite frequently hits API limits for Seqera container pulls; to resolve this, request a user token from Seqera as described in [troubleshooting.md](troubleshooting.md).

After tests finish, you may want to clean up by running `sudo rm -rf .nf-test/`.

> [!NOTE]
> Periodically delete docker images to free up space on your instance. Running the following command will delete all docker images on your system:
> ```
> docker kill $(docker ps -q) 2>/dev/null || true
> docker rm $(docker ps -a -q) 2>/dev/null || true
> docker rmi $(docker images -q) -f 2>/dev/null || true
> docker system prune -af --volumes
> ```

### Updating snapshots

Our test suite includes end-to-end tests of the `RUN` and `DOWNSTREAM` workflows that verify their outputs on test datasets have not changed unexpectedly. These tests use `nf-test`'s [snapshots](https://www.nf-test.com/docs/assertions/snapshots/).

The process of checking and updating snapshots when output has changed is a bit fiddly:
- These tests will fail when any output has changed, with the error message like:
```
Test [7677da69] 'RUN workflow output should match snapshot' 
  java.lang.RuntimeException: Different Snapshot:
  <md5 sums of previous output and new output>
```  
- First, make sure the changes are expected/desired:
    - Look at the md5 checksums to determine which files have changes; make sure they are what you expect.
    - Then, find the new output files and compare them to previous output files. Make sure the changes are expected based on your code changes.
        - Previous output files are in:
            - `test-data/gold-standard-results` (for short-read `RUN` test with the tag `run_output`)
            - `test-data/gold-standard-results-downstream` (for short-read `DOWNSTREAM` test with the tag `downstream`) 
            - `test-data/ont-ww-test-results` (for ONT `RUN` test with the tag `run_output_ont`)
        - New output files are in `.nf-test/tests/<hash>/output`. (`<hash>` is shown when the test runs; e.g., in the example error message above, the test hash begins with `7677da69`).
    - Once you are happy with the changes to the output:
        - Update output files in `test-data` by copying changed output files from the `.nf-test` directory to the appropriate location in `test-data`, uncompressing the files, and committing the changes.
        - Update `nf-test` snapshots by running `nf-test test <path to test that failed> --update-snapshot`; this will update the appropriate `*.snapshot` file in `tests/workflows`. Commit the changed snapshot file.


## GitHub issues
We use [GitHub issues](https://github.com/naobservatory/mgs-workflow/issues) to track any issues with the pipeline: bugs, cleanup tasks, and desired new features. 
Opening issues:
- All issues should be self-contained: the description should include enough detail to understand what needs to be done.
- We use labels to prioritize and track:
    - `enhancement`, `time&cost`, `bug`, and `documentation` to describe the type of issue.
    - `priority_1` (highest), `priority_2`, and `priority_3` to mark importance.
    - `in-progress` when actively being worked on, and `done` when resolved but changes are not yet merged to `master`.
         - only close issues when the fix/enhancement is merged to `master`.

## Pull requests (PRs)

### Creating branches and making changes
To contribute to the pipeline, start by creating a new branch off of `dev`. The branch should start with your name, followed by a short description of the feature you're working on, with hyphens (e.g. `will-post-hoc-integrate`).

>[!CAUTION]
> Do not make pull requests from the `dev` or `main` branches.

Please keep PRs small--a single feature or bugfix. For complex features, split across multiple PRs in small, logical chunks if at all possible. It's OK if a PR lays groundwork for a new feature but doesn't implement it yet.  

Feel free to use AI tools (Cursor, GitHub Copilot, Claude Code, etc.) to generate code and tests. The author is responsible for reviewing and understanding all AI-generated code before sending it for review.

### Sending PRs for review

1. **Write new tests** for the changes that you make using `nf-test` if those tests don't already exist. At the very least, these tests should check that the new implementation runs to completion; tests that also verify the output on the test dataset are strongly encouraged.
2. **Run all relevant tests locally** and make sure they pass. "Relevant" means: 1) Any tests of any process or workflow modified by the PR; 2) Any tests for any workflows that source any such process or workflow, and 3) Any tests that use any such process or workflow in setup.
    - You may run all existing tests as described in the "Testing" section above.
    - Or, you may identify relevant tests by recursively grepping for the process name in the `workflows`, `subworkflows`, and `tests` directories.
    - **Note which tests were run in your PR description.**
    - If you make any changes that affect the output of the pipeline, list/describe the changes that occurred in the pull request. 
3. **Update the `CHANGELOG.md` file** with the changes that you are making, and update the `pipeline-version.txt` file with the new version number.
    - More information on how to update the `CHANGELOG.md` file can be found [here](./versioning.md). Note that, before merging to `master`, version numbers should have the `-dev` suffix.
4. **Pass automated tests on GitHub Actions**. These run automatically when you open a pull request.
5. **Write a meaningful description** of your changes in the PR description and give it a meaningful title. 
    - In comments, feel free to flag any open questions or places where you need careful review. 
6. **Request review** from a maintainer on your changes. Current maintainers are jeffkaufman, willbradshaw, katherine-stansifer, and harmonbhasin. 
    - Make sure to assign the PR to the desired reviewer so that they see your PR (put them in the "Assignees" section on GitHub as well as in the "Reviewers" section).
7. To merge, you must **have an approving review** on your final changes, and all conversations must be resolved. After merging, please delete your branch! 


### Reviewing PRs

Reviewers should carefully read through all changes to code and tests, and make sure they understand them. Reviewers are not required to run/test code themselves.

Reviewers should make sure:
- The PR solves the desired issue in a logical and performant way.
- The code is clear, well-structured, and maintainable. 
- Tests and documentation have been updated appropriately. 
- The PR is a reasonable size and represents a coherent set of changes.
- The "before sending a PR" steps above have been followed.
- The guidelines in "Coding style guide" have been followed.

Make it clear whether your comments are:
- Items that must be addressed before submission
- Clarifying questions (you want a response to your comment but code changes are not required)
- Suggestions (author can take it or leave it)
- Issues that should be addressed, but could be in a later PR
    - Reviewer or author should create a GitHub issue to track anything in this category.

If you can't understand the code or the PR is overwhelmingly large, feel empowered to:
- Ask for the PR to be split into multiple smaller PRs. 
- Request more detailed comments and documentation.
- Make comments on the PR with clarifying questions ("Why do you need this function? What does this process do?")
- Request a meeting to walk through the PR in person.

When you are done with your review, assign the PR back to the author. 

As we require a review of final changes before merging, expect multiple rounds of review (with later rounds generally being very quick). After your final approving review, feel free to merge it and delete the branch!

## New releases

By default, all changes are made on individual branches, and merged into `dev`. Periodically, a collection of `dev` changes are merged to `master` as a new release. New releases are fairly frequent (historically, we have made a new release every 2-4 weeks).

Only a pipeline maintainer/member of the Nucleic Acid Observatory should go through the process below to author a new release. 

Process for a new release:
1. Merge `dev` into `staging`. (This does not require a PR or review.)
2. Update `CHANGELOG.md` and `pipeline-version.txt`.
    - See [versioning.md](./versioning.md) for information on our schematic versioning system.
    - Version numbers should increment by one with every new release (if `dev` has changes since the last release broken out into multiple smaller version bumps, combine them).
    - Remove any "-dev" suffix from version numbers. 
3. Double check that all documentation and tests are appropriately updated. 
4. Run *all* tests locally (see "Testing" section above).
5. Run the NAO's private test dataset through the `RUN` and `DOWNSTREAM` workflows.
    - Dataset can be found at `s3://nao-private-testing/300M-composite-dataset/`.
    - Check that it runs without errors; no need for additional validation of the output.
6. If the `INDEX` workflow is affected by the update, generate a new index from staging via the `INDEX` workflow and store it in `s3://nao-mgs-index` (private NAO bucket).
7. Make a PR to merge `staging` into `master`. Send it to a maintainer for review.
    - All changes should already have been reviewed as smaller PRs to `dev`; this review is a skim to make sure there are no unexpected changes, plus a careful review of the changelog and versioning. 
8. If any of the above surface any issues, fix them in `dev` and go back to step 1.
9. Once the PR to `master` is merged, create a new release (through the GitHub UI).
    - Title the release with the new version number.
    - For the release description, use the changelog entry for that version.
    - Tag the new release with the new version number (without any prepending "v").
10. Once the PR to `master` is merged, merge `staging` back into `dev`.








