#!/bin/bash

# Script to identify changed files compared to dev branch and run relevant tests
# For mgs-workflow repository based on new testing policy using nf-test

set -e

# Colors for better readability
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== MGS-Workflow Test Runner ===${NC}"
echo -e "${BLUE}Checking for changes compared to dev branch...${NC}"

# Get the list of changed files compared to dev branch
CHANGED_FILES=$(git diff --name-only origin/dev)

if [ -z "$CHANGED_FILES" ]; then
    echo -e "${YELLOW}No changes detected compared to dev branch.${NC}"
    exit 0
fi

echo -e "${GREEN}Changed files:${NC}"
echo "$CHANGED_FILES"
echo ""

# Arrays to store the affected components
declare -a AFFECTED_PROCESSES=()
declare -a AFFECTED_WORKFLOWS=()
declare -a AFFECTED_SUBWORKFLOWS=()
declare -a TESTS_TO_RUN=()

# Function to check if a string is in an array
contains() {
    local element="$1"
    shift
    for e in "$@"; do
        if [[ "$e" == "$element" ]]; then
            return 0
        fi
    done
    return 1
}

# Identify affected processes and workflows
for FILE in $CHANGED_FILES; do
    # Skip hidden files and directories
    if [[ "$FILE" == .* ]]; then
        continue
    fi

    # Check if file is a process
    if [[ "$FILE" == "modules/local/"*"/main.nf" ]]; then
        PROCESS_NAME=$(echo "$FILE" | sed 's/modules\/local\///g' | sed 's/\/main.nf//g')
        echo -e "${BLUE}Detected change in process: $PROCESS_NAME${NC}"
        if ! contains "$PROCESS_NAME" "${AFFECTED_PROCESSES[@]}"; then
            AFFECTED_PROCESSES+=("$PROCESS_NAME")
        fi

    # Check if file is a workflow
    elif [[ "$FILE" == "workflows/"*"/main.nf" ]]; then
        WORKFLOW_NAME=$(echo "$FILE" | sed 's/workflows\///g' | sed 's/\/main.nf//g')
        if ! contains "$WORKFLOW_NAME" "${AFFECTED_WORKFLOWS[@]}"; then
            AFFECTED_WORKFLOWS+=("$WORKFLOW_NAME")
        fi

    # Check if file is a subworkflow
    elif [[ "$FILE" == "subworkflows/local/"*"/main.nf" ]]; then
        SUBWORKFLOW_NAME=$(echo "$FILE" | sed 's/subworkflows\/local\///g' | sed 's/\/main.nf//g')
        if ! contains "$SUBWORKFLOW_NAME" "${AFFECTED_SUBWORKFLOWS[@]}"; then
            AFFECTED_SUBWORKFLOWS+=("$SUBWORKFLOW_NAME")
        fi
    fi
done

# Print identified components
if [ ${#AFFECTED_PROCESSES[@]} -gt 0 ]; then
    echo -e "${GREEN}Affected processes:${NC}"
    printf "  - %s\n" "${AFFECTED_PROCESSES[@]}"
fi

if [ ${#AFFECTED_WORKFLOWS[@]} -gt 0 ]; then
    echo -e "${GREEN}Affected workflows:${NC}"
    printf "  - %s\n" "${AFFECTED_WORKFLOWS[@]}"
fi

if [ ${#AFFECTED_SUBWORKFLOWS[@]} -gt 0 ]; then
    echo -e "${GREEN}Affected subworkflows:${NC}"
    printf "  - %s\n" "${AFFECTED_SUBWORKFLOWS[@]}"
fi

echo ""

# Find tests for directly modified components (category 1)
for PROCESS in "${AFFECTED_PROCESSES[@]}"; do
    # Find direct process tests - look for both *.nf.test and *_test.nf patterns
    echo -e "${BLUE}Looking for tests for process: $PROCESS${NC}"
    PROCESS_TESTS=$(find tests/modules/local -path "*${PROCESS}*/main.nf.test" -o -path "*${PROCESS}*/*_test.nf" 2>/dev/null || true)

    if [ -n "$PROCESS_TESTS" ]; then
        echo -e "${GREEN}Found tests:${NC}"
        echo "$PROCESS_TESTS"
        for TEST in $PROCESS_TESTS; do
            if ! contains "$TEST" "${TESTS_TO_RUN[@]}"; then
                TESTS_TO_RUN+=("$TEST")
            fi
        done
    else
        echo -e "${YELLOW}Warning: No direct tests found for process '$PROCESS'. You may need to create tests.${NC}"
        # Try with uppercase process name as fallback
        UPPERCASE_PROCESS=$(echo "$PROCESS" | tr '[:lower:]' '[:upper:]')
        PROCESS_TESTS=$(find tests/modules/local -path "*${UPPERCASE_PROCESS}*/main.nf.test" -o -path "*${UPPERCASE_PROCESS}*/*_test.nf" 2>/dev/null || true)

        if [ -n "$PROCESS_TESTS" ]; then
            echo -e "${GREEN}Found tests using uppercase process name:${NC}"
            echo "$PROCESS_TESTS"
            for TEST in $PROCESS_TESTS; do
                if ! contains "$TEST" "${TESTS_TO_RUN[@]}"; then
                    TESTS_TO_RUN+=("$TEST")
                fi
            done
        fi
    fi
done

for WORKFLOW in "${AFFECTED_WORKFLOWS[@]}"; do
    # Find direct workflow tests - look for both *.nf.test and *_test.nf patterns
    WORKFLOW_TESTS=$(find tests/workflows -path "*${WORKFLOW}*/main.nf.test" -o -path "*${WORKFLOW}*/*_test.nf" 2>/dev/null || true)
    if [ -n "$WORKFLOW_TESTS" ]; then
        for TEST in $WORKFLOW_TESTS; do
            if ! contains "$TEST" "${TESTS_TO_RUN[@]}"; then
                TESTS_TO_RUN+=("$TEST")
            fi
        done
    else
        echo -e "${YELLOW}Warning: No direct tests found for workflow '$WORKFLOW'. You may need to create tests.${NC}"
    fi
done

for SUBWORKFLOW in "${AFFECTED_SUBWORKFLOWS[@]}"; do
    # Find direct subworkflow tests - look for both *.nf.test and *_test.nf patterns
    SUBWORKFLOW_TESTS=$(find tests/subworkflows/local -path "*${SUBWORKFLOW}*/main.nf.test" -o -path "*${SUBWORKFLOW}*/*_test.nf" 2>/dev/null || true)
    if [ -n "$SUBWORKFLOW_TESTS" ]; then
        for TEST in $SUBWORKFLOW_TESTS; do
            if ! contains "$TEST" "${TESTS_TO_RUN[@]}"; then
                TESTS_TO_RUN+=("$TEST")
            fi
        done
    else
        echo -e "${YELLOW}Warning: No direct tests found for subworkflow '$SUBWORKFLOW'. You may need to create tests.${NC}"
    fi
done

# Find tests for workflows that source affected processes or subworkflows (category 2 & 3)
# This requires grep to find dependencies
for PROCESS in "${AFFECTED_PROCESSES[@]}"; do
    # Find workflows and subworkflows that use this process
    # Try both casing variants of the process name
    MIXED_CASE_PROCESS="$PROCESS"
    UPPERCASE_PROCESS=$(echo "$PROCESS" | tr '[:lower:]' '[:upper:]')

    DEPENDENT_WORKFLOWS=$(grep -r --include="*.nf" -e "${MIXED_CASE_PROCESS}" -e "${UPPERCASE_PROCESS}" workflows/ subworkflows/local/ 2>/dev/null || true)

    if [ -n "$DEPENDENT_WORKFLOWS" ]; then
        echo -e "${BLUE}Found workflows/subworkflows that depend on process '$PROCESS':${NC}"

        # Fix: Use a here-string instead of pipe to avoid subshell
        while read -r LINE; do
            FILE_PATH=$(echo "$LINE" | cut -d':' -f1)
            if [[ "$FILE_PATH" == "workflows/"*"/main.nf" ]]; then
                WORKFLOW_NAME=$(echo "$FILE_PATH" | sed 's/workflows\///g' | sed 's/\/main.nf//g')
                echo "  - Workflow: $WORKFLOW_NAME"

                # Find tests for this workflow
                WORKFLOW_TESTS=$(find tests/workflows -path "*${WORKFLOW_NAME}*/main.nf.test" -o -path "*${WORKFLOW_NAME}*/*_test.nf" 2>/dev/null || true)
                if [ -n "$WORKFLOW_TESTS" ]; then
                    for TEST in $WORKFLOW_TESTS; do
                        if ! contains "$TEST" "${TESTS_TO_RUN[@]}"; then
                            TESTS_TO_RUN+=("$TEST")
                        fi
                    done
                fi
            elif [[ "$FILE_PATH" == "subworkflows/local/"*"/main.nf" ]]; then
                SUBWORKFLOW_NAME=$(echo "$FILE_PATH" | sed 's/subworkflows\/local\///g' | sed 's/\/main.nf//g')
                echo "  - Subworkflow: $SUBWORKFLOW_NAME"

                # Find tests for this subworkflow
                SUBWORKFLOW_TESTS=$(find tests/subworkflows/local -path "*${SUBWORKFLOW_NAME}*/main.nf.test" -o -path "*${SUBWORKFLOW_NAME}*/*_test.nf" 2>/dev/null || true)
                if [ -n "$SUBWORKFLOW_TESTS" ]; then
                    for TEST in $SUBWORKFLOW_TESTS; do
                        if ! contains "$TEST" "${TESTS_TO_RUN[@]}"; then
                            TESTS_TO_RUN+=("$TEST")
                        fi
                    done
                fi
            fi
        done <<< "$DEPENDENT_WORKFLOWS"
    fi
done

for SUBWORKFLOW in "${AFFECTED_SUBWORKFLOWS[@]}"; do
    # Find workflows that use this subworkflow
    DEPENDENT_WORKFLOWS=$(grep -r --include="*.nf" "${SUBWORKFLOW}" workflows/ 2>/dev/null || true)

    if [ -n "$DEPENDENT_WORKFLOWS" ]; then
        echo -e "${BLUE}Found workflows that depend on subworkflow '$SUBWORKFLOW':${NC}"

        # Fix: Use a here-string instead of pipe to avoid subshell
        while read -r LINE; do
            FILE_PATH=$(echo "$LINE" | cut -d':' -f1)
            if [[ "$FILE_PATH" == "workflows/"*"/main.nf" ]]; then
                WORKFLOW_NAME=$(echo "$FILE_PATH" | sed 's/workflows\///g' | sed 's/\/main.nf//g')
                echo "  - Workflow: $WORKFLOW_NAME"

                # Find tests for this workflow
                WORKFLOW_TESTS=$(find tests/workflows -path "*${WORKFLOW_NAME}*/main.nf.test" -o -path "*${WORKFLOW_NAME}*/*_test.nf" 2>/dev/null || true)
                if [ -n "$WORKFLOW_TESTS" ]; then
                    for TEST in $WORKFLOW_TESTS; do
                        if ! contains "$TEST" "${TESTS_TO_RUN[@]}"; then
                            TESTS_TO_RUN+=("$TEST")
                        fi
                    done
                fi
            fi
        done <<< "$DEPENDENT_WORKFLOWS"
    fi
done

# Find tests that use affected components in setup (category 3)
for PROCESS in "${AFFECTED_PROCESSES[@]}"; do
    # Try both mixed case and uppercase process name
    MIXED_CASE_PROCESS="$PROCESS"
    UPPERCASE_PROCESS=$(echo "$PROCESS" | tr '[:lower:]' '[:upper:]')

    SETUP_TESTS=$(grep -r --include="*.nf.test" --include="*_test.nf" -e "${MIXED_CASE_PROCESS}" -e "${UPPERCASE_PROCESS}" tests/ 2>/dev/null || true)

    if [ -n "$SETUP_TESTS" ]; then
        echo -e "${BLUE}Found tests that use process '$PROCESS':${NC}"

        # Fix: Use a here-string instead of pipe to avoid subshell
        while read -r LINE; do
            FILE_PATH=$(echo "$LINE" | cut -d':' -f1)

            if [[ -f "$FILE_PATH" ]]; then
                echo "  - Test: $FILE_PATH"

                if ! contains "$FILE_PATH" "${TESTS_TO_RUN[@]}"; then
                    TESTS_TO_RUN+=("$FILE_PATH")
                fi
            fi
        done <<< "$SETUP_TESTS"
    fi
done

for SUBWORKFLOW in "${AFFECTED_SUBWORKFLOWS[@]}"; do
    SETUP_TESTS=$(grep -r --include="*.nf.test" --include="*_test.nf" "${SUBWORKFLOW}" tests/ 2>/dev/null || true)

    if [ -n "$SETUP_TESTS" ]; then
        echo -e "${BLUE}Found tests that use subworkflow '$SUBWORKFLOW' in setup:${NC}"
        while IFS= read -r LINE; do
            FILE_PATH=$(echo "$LINE" | cut -d':' -f1)

            if [[ -f "$FILE_PATH" ]]; then
                echo "  - Test: $FILE_PATH"

                if ! contains "$FILE_PATH" "${TESTS_TO_RUN[@]}"; then
                    TESTS_TO_RUN+=("$FILE_PATH")
                fi
            fi
        done <<< "$SETUP_TESTS"
    fi
done

echo ""

# Check if there are many changes that would warrant running the full test suite
CHANGE_COUNT=$(echo "$CHANGED_FILES" | wc -l)
if [ $CHANGE_COUNT -gt 10 ]; then
    echo -e "${YELLOW}This PR has a large number of changes ($CHANGE_COUNT files). Consider running the full test suite.${NC}"
    read -p "Do you want to run the full test suite instead? (y/n): " RUN_FULL
    if [[ "$RUN_FULL" == "y" || "$RUN_FULL" == "Y" ]]; then
        echo -e "${BLUE}Running full test suite...${NC}"
        nf-test test
        exit $?
    fi
fi

# Run identified tests
if [ ${#TESTS_TO_RUN[@]} -gt 0 ]; then
    echo -e "${GREEN}Tests to run:${NC}"
    printf "  - %s\n" "${TESTS_TO_RUN[@]}"

    echo ""
    echo -e "${BLUE}Running tests...${NC}"

    TEST_SUCCESS=true
    TEST_RESULTS=()

    for TEST in "${TESTS_TO_RUN[@]}"; do
        echo -e "${YELLOW}Running test: $TEST${NC}"

        # Run the test with nf-test
        if nf-test test $TEST; then
            TEST_RESULTS+=("✓ PASS: $TEST")
        else
            TEST_RESULTS+=("✗ FAIL: $TEST")
            TEST_SUCCESS=false
        fi
    done

    echo ""
    echo -e "${BLUE}Test Results:${NC}"
    printf "%s\n" "${TEST_RESULTS[@]}"

    if $TEST_SUCCESS; then
        echo -e "${GREEN}All tests passed successfully!${NC}"

        echo ""
        echo -e "${BLUE}PR Testing Summary:${NC}"
        echo -e "${GREEN}✓ Tested all modified processes and workflows${NC}"
        echo -e "${GREEN}✓ Tested all workflows that source modified processes/workflows${NC}"
        echo -e "${GREEN}✓ Tested all tests that use modified components in setup${NC}"
        echo -e "${GREEN}✓ All tests passed${NC}"

        echo ""
        echo -e "${YELLOW}Please include this summary in your PR description:${NC}"
        echo "## Testing"
        echo "I have run the following tests locally and they all pass:"
        for TEST in "${TESTS_TO_RUN[@]}"; do
            echo "- $TEST"
        done
        exit 0
    else
        echo -e "${RED}Some tests failed. Please fix the issues before submitting your PR.${NC}"
        exit 1
    fi
else
    echo -e "${YELLOW}No tests identified for the changed files. Please ensure you have appropriate test coverage.${NC}"
    exit 0
fi