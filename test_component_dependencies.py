#! /usr/bin/env python3
import argparse
import subprocess
import os

def find_component_dependencies(directory, component):

    directory = directory.rstrip("/")

    result = subprocess.run(
        ["grep", "-r", f"{component}", directory], capture_output=True, text=True
    )

    file_paths = set()
    for line in result.stdout.splitlines():
        if ":" in line:
            file_path = line.split(":", 1)[0]
            file_paths.add(file_path)

    return file_paths

def check_workflow_dependencies(workflow, script):

    workflow_path = f"workflows/{workflow}"
    script = script.replace("/main.nf", "")

    with open(workflow_path, "r") as f:
        workflow_content = f.read()

    if script in workflow_content:
        return True
    else:
        return False


def find_script_dependencies(directory, script):
    directory = directory.rstrip("/")
    script = script.replace("/main.nf", "")

    result = subprocess.run(
        ["grep", "-r", f"{script}", directory], capture_output=True, text=True
    )

    file_paths = set()
    for line in result.stdout.splitlines():
        if ":" in line:
            file_path = line.split(":", 1)[0]
            file_paths.add(file_path)

    return file_paths


def get_workflow_test(workflow):
    workflow_path = f"workflows/{workflow}"
    if "run_dev_se" in workflow_path:
        return "tests/workflows/run_dev.nf.test"
    else:
        base_name = os.path.basename(workflow_path)
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

    direct_tests = find_component_dependencies("tests/", component)
    tests_to_execute.update(direct_tests)

    dependent_subworkflows = set()
    for directory in os.listdir("subworkflows/"):
        dependent_subworkflows.update(find_component_dependencies(directory, component))


    subworkflow_tests = set()
    for dependent_subworkflow in dependent_subworkflows:
        found_tests = find_script_dependencies("tests/", dependent_subworkflow)
        subworkflow_tests.update(found_tests)

    workflow_tests = set()
    for workflow in os.listdir("workflows/"):
        for dependent_subworkflow in dependent_subworkflows:
            if check_workflow_dependencies(workflow, dependent_subworkflow):
                workflow_tests.add(get_workflow_test(workflow))


    tests_to_execute.update(subworkflow_tests)
    tests_to_execute.update(workflow_tests)

    for test in tests_to_execute:
        print(test)


    print("\n" + "=" * 80)
    print(f"Found {len(tests_to_execute)} total tests to run:")
    print("=" * 80)
    for test in sorted(tests_to_execute):
        print(f"   • {test}")

    print("\n" + "=" * 80)
    print("Running tests...")
    print("=" * 80)
    test_results = {}
    for test in sorted(tests_to_execute):
        print(f"\nRunning {test}")
        print("-" * 80)
        args = ["nf-test", "test", test]
        if verbose:
            args.append("--verbose")
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

    print(
        "Remember to run appropriate workflow tests (nf-test test tests/workflows/...)\n"
    )


if __name__ == "__main__":
    main()
