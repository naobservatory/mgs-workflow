#!/usr/bin/env python
"""
Test suite for the mark_ancestor_infections_single function.

This file tests the core logic of viral host infection status propagation from
child taxa to parent taxa. The function being tested determines a parent taxon's 
infection status based on its children's statuses, implementing the inference 
rules described in docs/annotation.md.

How to read this file:
1. Start with the test classes at the bottom - they show the actual test cases
2. Each test class covers a specific category of conditional statements in the function
3. The helper functions at the top generate test case combinations
4. Comments track which state combinations remain untested after each test

Key concepts:
- Parent nodes can only have initial states of UNRESOLVED (-1) or INCONSISTENT (0)
  (MATCH nodes are pre-marked as "checked" and never processed by this function)
- Child nodes can have any of the six states: UNRESOLVED, MATCH, INCONSISTENT,
  CONSISTENT, UNCLEAR, or MAYBE_INCONSISTENT
- The "Untested before/after" comments help track test coverage progression

Note: Additional integration testing exists in tests/modules/local/annotateVirusInfection/main.nf.test,
but comprehensive unit testing is being migrated to this file.

State meanings:
- MATCH (1): Confirmed to infect the host
- INCONSISTENT (0): Confirmed NOT to infect the host
- CONSISTENT (3): Likely infects (has MATCH descendants, no INCONSISTENT ones)
- UNCLEAR (2): Uncertain (has both MATCH and INCONSISTENT descendants)
- UNRESOLVED (-1): No data available (temporary state)
- MAYBE_INCONSISTENT (-2): Has INCONSISTENT descendants (temporary state)
"""

# =======================================================================
# Preamble
# =======================================================================

import pytest
import pandas as pd
import itertools
from typing import Callable, TypeAlias, Protocol

from annotate_viral_hosts import (
    mark_ancestor_infections_single,
    INCONSISTENT,
    MATCH,
    UNCLEAR,
    CONSISTENT,
    UNRESOLVED,
    MAYBE_INCONSISTENT,
)

# =======================================================================
# Type definitions and constants
# =======================================================================

State: TypeAlias = int
StateTuple: TypeAlias = tuple[State, ...]


class DataFrameFactory(Protocol):
    """Protocol for the dataframe factory fixture. Mainly used for type hinting."""

    def __call__(
        self,
        parent_status: int,
        child_statuses: dict[str, int],
        checked_states: dict[str, bool] | None = None,
    ) -> pd.DataFrame: ...


ALL_STATES: list[State] = [
    INCONSISTENT,
    MATCH,
    UNCLEAR,
    CONSISTENT,
    UNRESOLVED,
    MAYBE_INCONSISTENT,
]


# =======================================================================
# Test case generation helpers
# =======================================================================
def generate_combinations_with_state(
    state: State, num_children: int = 3
) -> list[StateTuple]:
    """Generate all combinations of states that include at least one instance of the given state."""
    all_combinations = itertools.product(ALL_STATES, repeat=num_children)
    return [
        case
        for case in all_combinations
        if state in case and len(set(case)) > 1  # Exclude uniform combinations
    ]


def generate_conflict_combinations(num_children: int = 3) -> list[StateTuple]:
    """Generate combinations of states with both inconsistent/maybe_inconsistent and match/consistent states."""
    inconsistent_group = {INCONSISTENT, MAYBE_INCONSISTENT}
    matched_group = {MATCH, CONSISTENT}
    all_combinations = itertools.product(ALL_STATES, repeat=num_children)
    return [
        combo
        for combo in all_combinations
        if any(item in inconsistent_group for item in combo)
        and any(item in matched_group for item in combo)
    ]


def generate_mixed_state_cases(
    child_states: list[State],
    parent_states: list[State],
    expected_logic: Callable[[tuple, int], int],
) -> list[StateTuple]:
    """Generate test cases for mixed child states with expected parent outcomes."""
    all_child_combos = itertools.product(child_states, repeat=3)
    # Filter out uniform combinations (all children have same state)
    mixed_combos = [combo for combo in all_child_combos if len(set(combo)) > 1]
    test_cases = []
    for children in mixed_combos:
        for initial_parent in parent_states:
            expected_parent = expected_logic(children, initial_parent)
            test_cases.append((*children, initial_parent, expected_parent))
    return test_cases


def inconsistent_propagation_logic(_: StateTuple, initial_parent: State) -> State:
    # We ignore children states because for this test case, all children are 
    # UNRESOLVED, MAYBE_INCONSISTENT, or INCONSISTENT. In this scenario,
    # only the parent's initial state determines the outcome.
    if initial_parent == INCONSISTENT:
        return INCONSISTENT
    return MAYBE_INCONSISTENT


# Generate propagation test cases
inconsistent_propagation_cases = generate_mixed_state_cases(
    child_states=[UNRESOLVED, INCONSISTENT, MAYBE_INCONSISTENT],
    parent_states=[UNRESOLVED, INCONSISTENT],
    expected_logic=inconsistent_propagation_logic,
)

consistent_propagation_cases = generate_mixed_state_cases(
    child_states=[UNRESOLVED, CONSISTENT, MATCH],
    parent_states=[UNRESOLVED, INCONSISTENT],
    expected_logic=lambda children, parent: CONSISTENT,
)


def generate_test_id(params: tuple) -> str:
    """Function to format test case IDs."""
    if len(params) == 5:
        c1, c2, c3, initial_p, expected_p = params
        return f"children({c1},{c2},{c3})-init({initial_p})-expect({expected_p})"
    else:
        return str(params)


# =======================================================================
# Fixtures
# =======================================================================


@pytest.fixture
def virus_taxonomy_tree() -> dict[str, set[str]]:
    """Create a test taxonomy tree with parent-child relationships."""

    return {
        "parent1": {"child1", "child2", "child3"},
        "child1": set(),
        "child2": set(),
        "child3": set(),
    }


@pytest.fixture
def dataframe_factory() -> DataFrameFactory:
    """Factory fixture to create test DataFrames with specified states."""

    def create_dataframe(
        parent_status: int,
        child_statuses: dict[str, int],
        checked_states: dict[str, bool] | None = None,
    ) -> pd.DataFrame:
        if checked_states is None:
            checked_states = {}
        rows = [
            {
                "taxid": "parent1",
                "status": parent_status,
                "checked": checked_states.get("parent1", False),
            }
        ]
        for child_id, status in child_statuses.items():
            rows.append(
                {
                    "taxid": child_id,
                    "status": status,
                    "checked": checked_states.get(child_id, True),
                }
            )

        df = pd.DataFrame(rows)
        df.set_index("taxid", inplace=True)
        return df

    return create_dataframe


# =======================================================================
# Test Classes
# =======================================================================


class TestUniformChildStates:
    # Untested before: All combinations
    @pytest.mark.parametrize(
        "child_state,expected_parent",
        [
            (UNRESOLVED, UNRESOLVED),
            (MATCH, MATCH),
            (INCONSISTENT, INCONSISTENT),
            (CONSISTENT, CONSISTENT),
            (UNCLEAR, UNCLEAR),
            (MAYBE_INCONSISTENT, MAYBE_INCONSISTENT),
        ],
    )
    def test_uniform_child_states_propagate_to_parent(
        self,
        dataframe_factory: DataFrameFactory,
        virus_taxonomy_tree: dict[str, set[str]],
        child_state: State,
        expected_parent: State,
    ) -> None:
        """If all children have the same state, the parent also has that state."""
        # Arrange
        child_statuses: dict[str, State] = {
            "child1": child_state,
            "child2": child_state,
            "child3": child_state,
        }
        df = dataframe_factory(
            parent_status=UNRESOLVED,
            child_statuses=child_statuses,
        )
        # Act
        result_df = mark_ancestor_infections_single("parent1", df, virus_taxonomy_tree)
        # Assert
        assert result_df.loc["parent1", "checked"]
        assert result_df.loc["parent1", "status"] == expected_parent

    # Untested after: All combinations where children have different states


class TestUnclearStatePropagation:
    # Untested before: All mixed child state combinations
    @pytest.mark.parametrize(
        "child1,child2,child3", generate_combinations_with_state(UNCLEAR)
    )
    def test_any_unclear_child_makes_parent_unclear(
        self,
        dataframe_factory: DataFrameFactory,
        virus_taxonomy_tree: dict[str, set[str]],
        child1: State,
        child2: State,
        child3: State,
    ) -> None:
        """Any child with UNCLEAR state should make parent state UNCLEAR."""
        # Arrange
        child_statuses: dict[str, State] = {
            "child1": child1,
            "child2": child2,
            "child3": child3,
        }
        df = dataframe_factory(parent_status=UNRESOLVED, child_statuses=child_statuses)
        # Act
        result_df = mark_ancestor_infections_single("parent1", df, virus_taxonomy_tree)
        # Assert
        assert result_df.loc["parent1", "checked"]
        assert result_df.loc["parent1", "status"] == UNCLEAR

    # Untested after: Mixed states without any UNCLEAR children

    @pytest.mark.parametrize("child1,child2,child3", generate_conflict_combinations())
    def test_conflicting_states_make_parent_unclear(
        self,
        dataframe_factory: DataFrameFactory,
        virus_taxonomy_tree: dict[str, set[str]],
        child1: State,
        child2: State,
        child3: State,
    ) -> None:
        """If there are children with both inconsistent/maybe_inconsistent and match/consistent states, make the parent UNCLEAR."""
        # Arrange
        child_statuses: dict[str, State] = {
            "child1": child1,
            "child2": child2,
            "child3": child3,
        }
        df = dataframe_factory(parent_status=UNRESOLVED, child_statuses=child_statuses)
        # Act
        result_df = mark_ancestor_infections_single("parent1", df, virus_taxonomy_tree)
        # Assert
        assert result_df.loc["parent1", "checked"]
        assert result_df.loc["parent1", "status"] == UNCLEAR

    # Untested after: Mixed states without UNCLEAR children AND without conflicting positive/negative evidence


class TestMixedStatePropagation:
    """Test propagation rules for mixed child states."""

    # Untested before: Mixed states without UNCLEAR children AND without conflicting positive/negative evidence

    @pytest.mark.parametrize(
        "child1,child2,child3,initial_parent,expected_parent",
        inconsistent_propagation_cases,
        ids=[generate_test_id(case) for case in inconsistent_propagation_cases],
    )
    def test_inconsistent_state_propagation(
        self,
        dataframe_factory: DataFrameFactory,
        virus_taxonomy_tree: dict[str, set[str]],
        child1: State,
        child2: State,
        child3: State,
        initial_parent: State,
        expected_parent: State,
    ) -> None:
        """Test propagation when children have inconsistent/maybe_inconsistent/unresolved states."""
        # Arrange
        child_statuses: dict[str, State] = {
            "child1": child1,
            "child2": child2,
            "child3": child3,
        }
        df = dataframe_factory(
            parent_status=initial_parent, child_statuses=child_statuses
        )
        # Act
        result_df = mark_ancestor_infections_single("parent1", df, virus_taxonomy_tree)
        # Assert
        assert result_df.loc["parent1", "checked"]
        assert result_df.loc["parent1", "status"] == expected_parent, (
            f"Expected parent status {expected_parent}, "
            f"but got {result_df.loc['parent1', 'status']}"
        )

    # Untested after: Mixed states with only positive evidence (MATCH/CONSISTENT/UNRESOLVED)

    @pytest.mark.parametrize(
        "child1,child2,child3,initial_parent,expected_parent",
        consistent_propagation_cases,
        ids=[generate_test_id(case) for case in consistent_propagation_cases],
    )
    def test_consistent_state_propagation(
        self,
        dataframe_factory: DataFrameFactory,
        virus_taxonomy_tree: dict[str, set[str]],
        child1: State,
        child2: State,
        child3: State,
        initial_parent: State,
        expected_parent: State,
    ) -> None:
        """Test propagation when children have match/consistent/unresolved states."""
        # Arrange
        child_statuses: dict[str, State] = {
            "child1": child1,
            "child2": child2,
            "child3": child3,
        }
        df = dataframe_factory(
            parent_status=initial_parent, child_statuses=child_statuses
        )
        # Act
        result_df = mark_ancestor_infections_single("parent1", df, virus_taxonomy_tree)
        # Assert
        assert result_df.loc["parent1", "checked"]
        assert result_df.loc["parent1", "status"] == expected_parent

    # Untested after: All major logic branches tested for mark_ancestor_infections_single


# =======================================================================
# Main execution
# =======================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
