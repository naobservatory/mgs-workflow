#! /usr/bin/env python3
import argparse
import subprocess
import os
import re
import sys
def find_dependency(directory, component):

    directory = directory.rstrip("/")
    component = component.replace("/main.nf", "")

    # Use grep with word boundaries to preempt substring matches.
    grep_results = subprocess.run(
        ["grep", "-r", f"\\b{component}\\b", directory], capture_output=True, text=True
    )

    # Relevant line triggers
    relevant_triggers = ["process", "include", "workflow"]

    # Extract file paths from relevant lines
    file_paths = {
        line.split(":", 1)[0]
        for line in grep_results.stdout.splitlines()
        if ":" in line
    }
    file_paths = set()
    for line in grep_results.stdout.splitlines():
        if any(trigger in line for trigger in relevant_triggers):
                file = line.split(":", 1)[0]
                file_paths.add(file)

    return file_paths


def workflow_uses_subworkflow(workflow, subworkflow_path):
    workflow_path = f"workflows/{workflow}"
    subworkflow_dir = subworkflow_path.replace("/main.nf", "")

    with open(workflow_path, "r") as f:
        workflow_content = f.read()

    # Use word boundaries to preempt substring matches
    if re.search(r'\b' + re.escape(subworkflow_dir) + r'\b', workflow_content):
        return True
    else:
        return False


def subworkflow_uses_subworkflow(subworkflow, dependent_subworkflow):
    subworkflow_path = f"subworkflows/local/{subworkflow}/main.nf"
    dependent_subworkflow_dir = dependent_subworkflow.replace("/main.nf", "")

    with open(subworkflow_path, "r") as f:
        subworkflow_content = f.read()

    # Use word boundaries to preempt substring matches
    if re.search(r'\b' + re.escape(dependent_subworkflow_dir) + r'\b', subworkflow_content):
        return True
    else:
        return False


def get_workflow_test(workflow):
    workflow_path = f"workflows/{workflow}"
    if re.search(r'\b' + re.escape("run_dev_se") + r'\b', workflow_path):
        return "tests/workflows/run_dev.nf.test"
    else:
        base_name = os.path.basename(workflow_path).split(".")[0]
        return f"tests/workflows/{base_name}.nf.test"


def parse_args():
    """Parse arguments and run the test finder."""
    parser = argparse.ArgumentParser(
        description="Find and run tests for specified components"
    )
    parser.add_argument(
        "-c",
        "--component",
        required=True,
        help="Name of component (process or subworkflow) to test",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose nf-test output",
        default=False,
    )

    args = parser.parse_args()
    component = args.component
    verbose = args.verbose
    return component, verbose


def main():
    component, verbose = parse_args()

    # Validate component name - only allow alphanumerics and underscores
    if not re.match(r'^[a-zA-Z0-9_]+$', component):
        raise ValueError(f"Error: Invalid component name '{component}'. Component names must contain only alphanumeric characters and underscores.")

    tests_to_execute = set()

    # Identifying tests that directly use the component
    direct_tests = find_dependency("tests/", component)
    tests_to_execute.update(direct_tests)

    print("=" * 72)
    print(f"Found {len(direct_tests)} direct tests:")
    print("=" * 72)
    for test in direct_tests:
        print(f"   • {test}")
    print()

    # Find subworkflows that directly use the component
    dependent_subworkflows = set(find_dependency("subworkflows/", component))

    # Find transitive dependencies - subworkflows that use other dependent subworkflows
    transitive_dependencies = set()
    for subworkflow in os.listdir("subworkflows/local/"):
        for dependent_subworkflow in dependent_subworkflows:
            if subworkflow_uses_subworkflow(subworkflow, dependent_subworkflow):
                transitive_dependencies.add(subworkflow)

    # Combine all dependent subworkflows
    dependent_subworkflows.update(transitive_dependencies)

    # Identifying tests that use the affected subworkflows
    subworkflow_tests = set()
    for dependent_subworkflow in dependent_subworkflows:
        tests = find_dependency("tests/subworkflows/", dependent_subworkflow)
        subworkflow_tests.update(tests)

    print("=" * 72)
    print(f"Found {len(subworkflow_tests)} tests that depend on affected subworkflows:")
    print("=" * 72)
    for test in subworkflow_tests:
        print(f"   • {test}")
    print()

    # Identifying workflows that use the component
    dependent_workflows = set()
    dependent_workflows.update(find_dependency("workflows/", component))

    # Identifying workflows that use affected subworkflows
    for workflow in os.listdir("workflows/"):
        for dependent_subworkflow in dependent_subworkflows:
            if workflow_uses_subworkflow(workflow, dependent_subworkflow):
                dependent_workflows.add(workflow)

    # Identifying tests that use the affected workflows
    workflow_tests = set()
    for workflow in dependent_workflows:
        workflow_tests.add(get_workflow_test(workflow))

    print("=" * 72)
    print(f"Found {len(workflow_tests)} dependent workflow tests:")
    print("=" * 72)
    for test in workflow_tests:
        print(f"   • {test}")
    print()

    tests_to_execute.update(subworkflow_tests)
    tests_to_execute.update(workflow_tests)

    print("=" * 72)
    print(f"Identified {len(tests_to_execute)} tests to execute. Running...")

    # Run all tests in a single nf-test command
    cmd = ["nf-test", "test"] + sorted(tests_to_execute)
    if verbose:
        cmd.append("--verbose")
        print(f"\nRunning test files with verbose output:")
    else:
        print(f"\nRunning test files:")

    for test in sorted(tests_to_execute):
        print(f"   • {test}")
    print("-" * 72)

    result = subprocess.run(cmd)


if __name__ == "__main__":
    main()
