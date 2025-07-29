#!/usr/bin/env python

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
    """Protocol for the dataframe factory fixture."""

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


def inconsistent_propagation_logic(
    children: StateTuple, initial_parent: State
) -> State:
    if initial_parent == INCONSISTENT and (
        INCONSISTENT in children or MAYBE_INCONSISTENT in children
    ):
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
            checked_states={"parent1": False},
        )
        # Act
        result_df = mark_ancestor_infections_single("parent1", df, virus_taxonomy_tree)
        # Assert
        assert result_df.loc["parent1", "checked"]
        assert result_df.loc["parent1", "status"] == expected_parent


class TestUnclearStatePropagation:
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


class TestMixedStatePropagation:
    """Test propagation rules for mixed child states."""

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


# =======================================================================
# Main execution
# =======================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
