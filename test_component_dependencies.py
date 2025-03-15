#! /usr/bin/env python3
import argparse
import subprocess


def find_files(directory, component):
    """
    Recursively grep files in a directory, and return file paths.
    """
    directory = directory.rstrip("/")
    component = component.replace("/main.nf", "")

    result = subprocess.run(
        ["grep", "-r", f"\\b{component}\\b", directory], capture_output=True, text=True
    )

    file_paths = set()
    for line in result.stdout.splitlines():
        if ":" in line:
            file_path = line.split(":", 1)[0]
            file_paths.add(file_path)

    return file_paths


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

    # Find tests that directly reference the component
    print("\n" + "=" * 80)
    print(f"Finding tests for component: {component}...")
    print("=" * 80)
    direct_tests = find_files("tests/", component)
    tests_to_execute.update(direct_tests)

    if direct_tests:
        print(f"\nFound {len(direct_tests)} tests that directly reference {component}:")
        for test in sorted(direct_tests):
            print(f"   • {test}")

    # Find subworkflows that use the component
    print("\n" + "=" * 80)
    print(f"Finding subworkflows that use {component}...")
    print("=" * 80)
    dependent_components = set()
    for directory in ["subworkflows/"]:
        dependent_components.update(find_files(directory, component))

    if dependent_components:
        print(f"\nFound {len(dependent_components)} subworkflows that use {component}:")
        for dependent_component in sorted(dependent_components):
            print(f"   • {dependent_component}")

    # Find tests for those subworkflows
    print("\n" + "=" * 80)
    print(f"Finding tests for subworkflows that use {component}...")
    print("=" * 80)
    indirect_tests = set()
    for dependent_component in dependent_components:
        found_tests = find_files("tests/", dependent_component)
        indirect_tests.update(found_tests)

    if indirect_tests:
        print(
            f"\nFound {len(indirect_tests)} tests for subworkflows that use {component}:"
        )
        for test in sorted(indirect_tests):
            print(f"   • {test}")

    # Combine all tests
    tests_to_execute.update(indirect_tests)

    # Display results
    print("\n" + "=" * 80)
    print(f"Found {len(tests_to_execute)} total tests to run:")
    print("=" * 80)
    for test in sorted(tests_to_execute):
        print(f"   • {test}")

    # Run the tests
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

    # Print summary with colored checkmarks
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
