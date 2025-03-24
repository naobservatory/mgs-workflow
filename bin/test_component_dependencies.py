#! /usr/bin/env python3
import argparse
import subprocess
import os


def find_dependency(directory, component):

    directory = directory.rstrip("/")
    component = component.replace("/main.nf", "")

    grep_results = subprocess.run(
        ["grep", "-r", component, directory], capture_output=True, text=True
    )

    file_paths = {
        line.split(":", 1)[0]
        for line in grep_results.stdout.splitlines()
        if ":" in line
    }

    return file_paths


def workflow_uses_subworkflow(workflow, subworkflow_path):

    workflow_path = f"workflows/{workflow}"
    subworkflow_dir = subworkflow_path.replace("/main.nf", "")

    with open(workflow_path, "r") as f:
        workflow_content = f.read()

    if subworkflow_dir in workflow_content:
        return True
    else:
        return False


def get_workflow_test(workflow):
    workflow_path = f"workflows/{workflow}"
    if "run_dev_se" in workflow_path:
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

    tests_to_execute = set()

    direct_tests = find_dependency("tests/", component)
    tests_to_execute.update(direct_tests)

    print("=" * 80)
    print(f"Found {len(direct_tests)} direct tests:")
    print("=" * 80)
    for test in direct_tests:
        print(f"   • {test}")
    print()
    # print("\n")

    dependent_subworkflows = set()
    subworkflow_tests = set()

    dependent_subworkflows.update(find_dependency("subworkflows/", component))
    for dependent_subworkflow in dependent_subworkflows:
        tests = find_dependency("tests/", dependent_subworkflow)
        subworkflow_tests.update(tests)

    print("=" * 80)
    print(f"Found {len(subworkflow_tests)} dependent subworkflow tests:")
    print("=" * 80)
    for test in subworkflow_tests:
        print(f"   • {test}")
    print()

    workflow_tests = set()
    for workflow in os.listdir("workflows/"):
        for dependent_subworkflow in dependent_subworkflows:
            if workflow_uses_subworkflow(workflow, dependent_subworkflow):
                workflow_tests.add(get_workflow_test(workflow))

    print("=" * 80)
    print(f"Found {len(workflow_tests)} dependent workflow tests:")
    print("=" * 80)
    for test in workflow_tests:
        print(f"   • {test}")
    print()

    tests_to_execute.update(subworkflow_tests)
    tests_to_execute.update(workflow_tests)

    print("=" * 80)
    print(f"Identified {len(tests_to_execute)} tests to execute. Running...")

    test_results = {}
    for test in sorted(tests_to_execute):
        args = ["nf-test", "test", test]
        if verbose:
            args.append("--verbose")
        print(f"\nRunning \"{' '.join(args)}\"")
        print("-" * 80)
        result = subprocess.run(args)

        test_results[test] = result.returncode == 0

    print("\n" + "=" * 80)
    print("Test Summary")
    print("=" * 80)
    for test, passed in test_results.items():
        if passed:
            status = "\033[92m✓\033[0m"  # Green checkmark
        else:
            status = "\033[91m✗\033[0m"  # Red X
        print(f"   {status} {test}")

    total = len(test_results)
    passed = sum(test_results.values())
    failed = total - passed

    print("\n" + "-" * 40)
    print(f"Total:  {total} tests")
    print(f"Passed: {passed} tests")
    if failed > 0:
        print(f"\033[91mFailed: {failed} tests\033[0m")
    else:
        print(f"Failed: {failed} tests")
    print("-" * 40 + "\n")


if __name__ == "__main__":
    main()
