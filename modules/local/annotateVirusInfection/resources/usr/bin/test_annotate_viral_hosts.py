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
    mark_direct_infections,
    expand_taxid,
    build_virus_tree,
    mark_descendant_infections,
    add_descendants,
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

# =======================================================================
# Tests for mark_ancestor_infections_single
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
# Tests for mark_direct_infections
# =======================================================================

class TestMarkDirectInfections:
    """Test the mark_direct_infections function."""

    def test_virus_has_matching_hosts(self) -> None:
        """Test that virus with matching hosts returns MATCH."""
        # Arrange
        virus_taxids = pd.Series(["virus1", "virus2", "virus3"])
        host_taxids = {"host1", "host2", "host3"}
        virus_host_mapping = {
            "virus1": {"host1", "host4"},  # Has matching host (host1)
            "virus2": {"host2", "host3"},  # Has matching hosts (host2, host3)
            "virus3": {"host4", "host5"},  # No matching hosts
        }
        
        # Act
        result = mark_direct_infections(virus_taxids, host_taxids, virus_host_mapping)
        
        # Assert
        assert result["virus1"] == MATCH
        assert result["virus2"] == MATCH
        assert result["virus3"] == INCONSISTENT

    def test_virus_has_non_matching_hosts(self) -> None:
        """Test that virus with non-matching hosts returns INCONSISTENT."""
        # Arrange
        virus_taxids = pd.Series(["virus1", "virus2"])
        host_taxids = {"host1", "host2"}
        virus_host_mapping = {
            "virus1": {"host3", "host4"},  # No matching hosts
            "virus2": {"host5"},           # No matching hosts
        }
        
        # Act
        result = mark_direct_infections(virus_taxids, host_taxids, virus_host_mapping)
        
        # Assert
        assert result["virus1"] == INCONSISTENT
        assert result["virus2"] == INCONSISTENT

    def test_virus_not_in_mapping(self) -> None:
        """Test that virus not in mapping returns UNRESOLVED."""
        # Arrange
        virus_taxids = pd.Series(["virus1", "virus2", "virus3"])
        host_taxids = {"host1", "host2"}
        virus_host_mapping = {
            "virus1": {"host1"},  # Has matching host
        }
        # virus2 and virus3 are not in mapping
        
        # Act
        result = mark_direct_infections(virus_taxids, host_taxids, virus_host_mapping)
        
        # Assert
        assert result["virus1"] == MATCH
        assert result["virus2"] == UNRESOLVED
        assert result["virus3"] == UNRESOLVED

    def test_empty_inputs(self) -> None:
        """Test that empty virus series is handled gracefully."""
        # Arrange
        virus_taxids = pd.Series([], dtype=str)
        host_taxids = {"host1", "host2"}
        virus_host_mapping = {"virus1": {"host1"}}
        
        # Act
        result = mark_direct_infections(virus_taxids, host_taxids, virus_host_mapping)
        
        # Assert
        assert len(result) == 0
        assert isinstance(result, pd.Series)

    def test_empty_host_taxids(self) -> None:
        """Test behavior with empty host taxids set."""
        # Arrange
        virus_taxids = pd.Series(["virus1", "virus2"])
        host_taxids = set()  # Empty set
        virus_host_mapping = {
            "virus1": {"host1", "host2"},
            "virus2": {"host3"},
        }
        
        # Act
        result = mark_direct_infections(virus_taxids, host_taxids, virus_host_mapping)
        
        # Assert
        # All viruses should be INCONSISTENT since no hosts match empty set
        assert result["virus1"] == INCONSISTENT
        assert result["virus2"] == INCONSISTENT


# =======================================================================
# Tests for expand_taxid
# =======================================================================

class TestExpandTaxid:
    """Test the expand_taxid function."""

    def test_linear_taxonomy(self) -> None:
        """Test expansion of a linear chain of descendants."""
        # Arrange - Create a linear taxonomy: parent -> child -> grandchild
        nodes_data = [
            {"taxid": "1", "parent_taxid": "0"},  # Root
            {"taxid": "2", "parent_taxid": "1"},  # Child of 1
            {"taxid": "3", "parent_taxid": "2"},  # Grandchild of 1
            {"taxid": "4", "parent_taxid": "3"},  # Great-grandchild of 1
        ]
        nodes = pd.DataFrame(nodes_data)
        
        # Act
        result = expand_taxid("1", nodes)
        
        # Assert
        expected = {"1", "2", "3", "4"}
        assert result == expected

    def test_branching_taxonomy(self) -> None:
        """Test expansion with multiple descendant branches."""
        # Arrange - Create a branching taxonomy
        nodes_data = [
            {"taxid": "1", "parent_taxid": "0"},  # Root
            {"taxid": "2", "parent_taxid": "1"},  # First child
            {"taxid": "3", "parent_taxid": "1"},  # Second child
            {"taxid": "4", "parent_taxid": "2"},  # Grandchild via first child
            {"taxid": "5", "parent_taxid": "3"},  # Grandchild via second child
            {"taxid": "6", "parent_taxid": "3"},  # Another grandchild via second child
        ]
        nodes = pd.DataFrame(nodes_data)
        
        # Act
        result = expand_taxid("1", nodes)
        
        # Assert
        expected = {"1", "2", "3", "4", "5", "6"}
        assert result == expected

    def test_leaf_node(self) -> None:
        """Test that leaf node returns only itself."""
        # Arrange
        nodes_data = [
            {"taxid": "1", "parent_taxid": "0"},
            {"taxid": "2", "parent_taxid": "1"},
            {"taxid": "3", "parent_taxid": "1"},
            {"taxid": "4", "parent_taxid": "2"},  # Leaf node
        ]
        nodes = pd.DataFrame(nodes_data)
        
        # Act
        result = expand_taxid("4", nodes)  # Query leaf node
        
        # Assert
        assert result == {"4"}

    def test_taxid_not_in_nodes(self) -> None:
        """Test that taxid not in nodes returns only itself."""
        # Arrange
        nodes_data = [
            {"taxid": "1", "parent_taxid": "0"},
            {"taxid": "2", "parent_taxid": "1"},
        ]
        nodes = pd.DataFrame(nodes_data)
        
        # Act
        result = expand_taxid("999", nodes)  # Non-existent taxid
        
        # Assert
        assert result == {"999"}

    def test_complex_tree(self) -> None:
        """Test expansion in a more complex tree structure."""
        # Arrange - Create a complex tree
        nodes_data = [
            {"taxid": "1", "parent_taxid": "0"},
            # First branch
            {"taxid": "10", "parent_taxid": "1"},
            {"taxid": "11", "parent_taxid": "10"},
            {"taxid": "12", "parent_taxid": "10"},
            # Second branch
            {"taxid": "20", "parent_taxid": "1"},
            {"taxid": "21", "parent_taxid": "20"},
            {"taxid": "22", "parent_taxid": "20"},
            # Deep branch from 21
            {"taxid": "210", "parent_taxid": "21"},
            {"taxid": "211", "parent_taxid": "210"},
        ]
        nodes = pd.DataFrame(nodes_data)
        
        # Act - Expand from different starting points
        result_root = expand_taxid("1", nodes)
        result_branch = expand_taxid("20", nodes)
        
        # Assert
        expected_root = {"1", "10", "11", "12", "20", "21", "22", "210", "211"}
        expected_branch = {"20", "21", "22", "210", "211"}
        assert result_root == expected_root
        assert result_branch == expected_branch


# =======================================================================
# Tests for build_virus_tree
# =======================================================================

class TestBuildVirusTree:
    """Test the build_virus_tree function."""

    def test_valid_input(self) -> None:
        """Test creating correct parent→children mapping from valid input."""
        # Arrange
        viral_taxa_data = [
            {"taxid": "1", "parent_taxid": "0"},  # Root has parent 0
            {"taxid": "2", "parent_taxid": "1"},  # Child of 1
            {"taxid": "3", "parent_taxid": "1"},  # Another child of 1
            {"taxid": "4", "parent_taxid": "2"},  # Child of 2
            {"taxid": "5", "parent_taxid": "2"},  # Another child of 2
        ]
        viral_taxa_df = pd.DataFrame(viral_taxa_data)
        
        # Act
        result = build_virus_tree(viral_taxa_df)
        
        # Assert
        expected = {
            "0": {"1"},        # Parent 0 has child 1
            "1": {"2", "3"},  # Parent 1 has children 2 and 3
            "2": {"4", "5"},  # Parent 2 has children 4 and 5
        }
        # Check all expected keys are present
        for parent, children in expected.items():
            assert parent in result
            assert result[parent] == children
        
        # Verify that nodes without children return empty sets
        assert result["3"] == set()
        assert result["4"] == set()
        assert result["5"] == set()

    def test_empty_dataframe(self) -> None:
        """Test that empty DataFrame returns empty defaultdict."""
        # Arrange
        viral_taxa_df = pd.DataFrame(columns=["taxid", "parent_taxid"])
        
        # Act
        result = build_virus_tree(viral_taxa_df)
        
        # Assert
        # Should be an empty defaultdict
        assert len(result) == 0
        # But still should behave like defaultdict
        assert result["any_key"] == set()

    def test_multiple_children(self) -> None:
        """Test parent with multiple children."""
        # Arrange
        viral_taxa_data = [
            {"taxid": "10", "parent_taxid": "1"},
            {"taxid": "20", "parent_taxid": "1"},
            {"taxid": "30", "parent_taxid": "1"},
            {"taxid": "40", "parent_taxid": "1"},
            {"taxid": "50", "parent_taxid": "1"},
        ]
        viral_taxa_df = pd.DataFrame(viral_taxa_data)
        
        # Act
        result = build_virus_tree(viral_taxa_df)
        
        # Assert
        assert result["1"] == {"10", "20", "30", "40", "50"}
        # All children should have empty sets
        for child_id in ["10", "20", "30", "40", "50"]:
            assert result[child_id] == set()

    def test_complex_hierarchy(self) -> None:
        """Test a more complex viral taxonomy hierarchy."""
        # Arrange - Create a complex tree structure
        viral_taxa_data = [
            # Top level
            {"taxid": "100", "parent_taxid": "root"},
            {"taxid": "200", "parent_taxid": "root"},
            # Second level
            {"taxid": "110", "parent_taxid": "100"},
            {"taxid": "120", "parent_taxid": "100"},
            {"taxid": "210", "parent_taxid": "200"},
            # Third level
            {"taxid": "111", "parent_taxid": "110"},
            {"taxid": "112", "parent_taxid": "110"},
            {"taxid": "211", "parent_taxid": "210"},
            {"taxid": "212", "parent_taxid": "210"},
            {"taxid": "213", "parent_taxid": "210"},
        ]
        viral_taxa_df = pd.DataFrame(viral_taxa_data)
        
        # Act
        result = build_virus_tree(viral_taxa_df)
        
        # Assert
        assert result["root"] == {"100", "200"}
        assert result["100"] == {"110", "120"}
        assert result["200"] == {"210"}
        assert result["110"] == {"111", "112"}
        assert result["210"] == {"211", "212", "213"}
        # Leaf nodes
        assert result["120"] == set()
        assert result["111"] == set()
        assert result["112"] == set()
        assert result["211"] == set()
        assert result["212"] == set()
        assert result["213"] == set()


# =======================================================================
# Tests for mark_descendant_infections
# =======================================================================

class TestMarkDescendantInfections:
    """Test mark_descendant_infections function.
    
    This function propagates infection statuses from parent nodes to their descendants.
    The propagation follows a specific priority order:
    1. MATCH propagation (overrides everything)
    2. INCONSISTENT propagation (overrides everything except MATCH)
    3. CONSISTENT propagation (doesn't override MATCH)
    4. UNCLEAR propagation (only affects UNRESOLVED and MAYBE_INCONSISTENT)
    
    After this function runs:
    - No nodes should have UNRESOLVED status
    - No nodes should have MAYBE_INCONSISTENT status
    
    Valid Parent-Child State Combinations
    Based on the logic in the code, these are the valid parent-child combinations that can exist when mark_descendant_infections is called:
    
    Parent: MATCH
    - Children can be: Any state
    - Note: MATCH parents are marked as "checked" early and skip mark_ancestor_infections_single
    
    Parent: INCONSISTENT
    - If all children are INCONSISTENT: Parent inherits INCONSISTENT (uniform rule)
    - Otherwise: Children can be INCONSISTENT, MAYBE_INCONSISTENT, or UNRESOLVED
    - Cannot have: MATCH or CONSISTENT children (would make parent UNCLEAR instead)
    
    Parent: UNCLEAR
    - Must have one of:
      - At least one UNCLEAR child, OR
      - Mix of (INCONSISTENT/MAYBE_INCONSISTENT) and (MATCH/CONSISTENT) children
    - Cannot have: All children with same state (would trigger uniform rule)
    
    Parent: CONSISTENT
    - Children must be some mix of: MATCH, CONSISTENT, and/or UNRESOLVED only
    - Cannot have: INCONSISTENT, MAYBE_INCONSISTENT, or UNCLEAR children
    - Cannot have: All children with same state (would trigger uniform rule)
    
    Parent: UNRESOLVED
    - If parent is UNRESOLVED after mark_ancestor_infections, children must all be UNRESOLVED
    - This follows from the logic: parent only stays UNRESOLVED if all children are UNRESOLVED
    
    Parent: MAYBE_INCONSISTENT
    - Must have: At least one INCONSISTENT or MAYBE_INCONSISTENT child
    - Parent's original state was UNRESOLVED (before mark_ancestor_infections_single)
    - Cannot have: Only MATCH/CONSISTENT children (would make parent CONSISTENT)
    """
    
    @pytest.fixture
    def simple_tree(self) -> dict[str, set[str]]:
        """Create a simple tree structure for testing.
        
        Structure:
        root
        ├── v1
        │   ├── v11
        │   └── v12
        └── v2
            └── v21
        """
        return {
            "root": {"v1", "v2"},
            "v1": {"v11", "v12"},
            "v2": {"v21"},
            "v11": set(),
            "v12": set(),
            "v21": set(),
        }
    
    @pytest.fixture
    def complex_tree(self) -> dict[str, set[str]]:
        """Create a more complex tree structure for testing.
        
        Structure:
        root
        ├── v1
        │   ├── v11
        │   │   └── v111
        │   └── v12
        ├── v2
        │   └── v21
        │       ├── v211
        │       └── v212
        └── v3
        """
        return {
            "root": {"v1", "v2", "v3"},
            "v1": {"v11", "v12"},
            "v2": {"v21"},
            "v3": set(),
            "v11": {"v111"},
            "v12": set(),
            "v21": {"v211", "v212"},
            "v111": set(),
            "v211": set(),
            "v212": set(),
        }

    def test_match_propagation_simple(self, simple_tree: dict[str, set[str]]) -> None:
        """Test that MATCH status propagates to all descendants.
        
        Initial state:
        root (UNCLEAR)
        ├── v1 (MATCH)
        │   ├── v11 (UNRESOLVED)
        │   └── v12 (INCONSISTENT)
        └── v2 (INCONSISTENT)
            └── v21 (UNRESOLVED)
        
        Expected after propagation:
        root (UNCLEAR)
        ├── v1 (MATCH)
        │   ├── v11 (MATCH)         # Inherited from MATCH parent
        │   └── v12 (MATCH)         # MATCH overrides INCONSISTENT
        └── v2 (INCONSISTENT)
            └── v21 (INCONSISTENT)  # Inherited from INCONSISTENT parent
        """
        # Arrange
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": MATCH,
            "v11": UNRESOLVED,
            "v12": INCONSISTENT,
            "v2": INCONSISTENT,
            "v21": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(simple_tree, statuses.copy())
        
        # Assert - v1 and all its descendants should be MATCH
        assert result["v1"] == MATCH
        assert result["v11"] == MATCH
        assert result["v12"] == MATCH
        # v2 and v21 should remain INCONSISTENT (overridden by INCONSISTENT propagation)
        assert result["v2"] == INCONSISTENT
        assert result["v21"] == INCONSISTENT
        # root should remain UNCLEAR
        assert result["root"] == UNCLEAR

    def test_match_propagation_overrides_all(self, simple_tree: dict[str, set[str]]) -> None:
        """Test that MATCH propagation overrides all other statuses.
        
        Initial state:
        root (MATCH)
        ├── v1 (INCONSISTENT)
        │   ├── v11 (UNCLEAR)
        │   └── v12 (CONSISTENT)
        └── v2 (MAYBE_INCONSISTENT)
            └── v21 (UNRESOLVED)
        
        Expected after propagation:
        All nodes become MATCH (MATCH propagates from root and overrides everything)
        """
        # Arrange - Set all descendants to different statuses
        statuses = pd.Series({
            "root": MATCH,
            "v1": INCONSISTENT,
            "v11": UNCLEAR,
            "v12": CONSISTENT,
            "v2": MAYBE_INCONSISTENT,
            "v21": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(simple_tree, statuses.copy())
        
        # Assert - Everything should be MATCH
        assert all(result == MATCH)

    def test_inconsistent_propagation(self, simple_tree: dict[str, set[str]]) -> None:
        """Test that INCONSISTENT status propagates to descendants.
        
        Initial state:
        root (UNCLEAR)
        ├── v1 (INCONSISTENT)
        │   ├── v11 (UNRESOLVED)
        │   └── v12 (MAYBE_INCONSISTENT)
        └── v2 (CONSISTENT)
            └── v21 (UNRESOLVED)
        
        Expected after propagation:
        root (UNCLEAR)
        ├── v1 (INCONSISTENT)
        │   ├── v11 (INCONSISTENT)  # Inherited from INCONSISTENT parent
        │   └── v12 (INCONSISTENT)  # INCONSISTENT overrides MAYBE_INCONSISTENT
        └── v2 (CONSISTENT)
            └── v21 (CONSISTENT)     # Inherited from CONSISTENT parent
        """
        # Arrange
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": INCONSISTENT,
            "v11": UNRESOLVED,
            "v12": MAYBE_INCONSISTENT,
            "v2": CONSISTENT,
            "v21": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(simple_tree, statuses.copy())
        
        # Assert
        assert result["v1"] == INCONSISTENT
        assert result["v11"] == INCONSISTENT
        assert result["v12"] == INCONSISTENT
        assert result["v2"] == CONSISTENT
        assert result["v21"] == CONSISTENT  # Inherits from CONSISTENT parent
        assert result["root"] == UNCLEAR

    def test_consistent_propagation(self, simple_tree: dict[str, set[str]]) -> None:
        """Test that CONSISTENT status propagates correctly.
        
        Initial state:
        root (CONSISTENT)
        ├── v1 (MATCH)
        │   ├── v11 (UNRESOLVED)
        │   └── v12 (UNRESOLVED)
        └── v2 (UNRESOLVED)
            └── v21 (UNRESOLVED)
        
        Expected after propagation:
        root (CONSISTENT)
        ├── v1 (MATCH)              # Preserved (not overridden)
        │   ├── v11 (MATCH)         # Inherits from MATCH parent (higher priority)
        │   └── v12 (MATCH)         # Inherits from MATCH parent (higher priority)
        └── v2 (CONSISTENT)         # Inherits from CONSISTENT root
            └── v21 (CONSISTENT)    # Inherits from CONSISTENT parent
        """
        # Arrange
        statuses = pd.Series({
            "root": CONSISTENT,
            "v1": MATCH,
            "v11": UNRESOLVED,
            "v12": UNRESOLVED,
            "v2": UNRESOLVED,
            "v21": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(simple_tree, statuses.copy())
        
        # Assert
        assert result["root"] == CONSISTENT
        assert result["v1"] == MATCH  # MATCH is not overridden
        assert result["v11"] == MATCH  # Inherits from MATCH parent first
        assert result["v12"] == MATCH  # Inherits from MATCH parent first
        assert result["v2"] == CONSISTENT  # Inherits from CONSISTENT root
        assert result["v21"] == CONSISTENT  # Inherits from CONSISTENT parent

    def test_consistent_propagation_preserves_match(self, complex_tree: dict[str, set[str]]) -> None:
        """Test that CONSISTENT propagation doesn't override MATCH status.
        
        Initial state:
        root (CONSISTENT)
        ├── v1 (CONSISTENT)
        │   ├── v11 (MATCH)
        │   │   └── v111 (UNRESOLVED)
        │   └── v12 (UNRESOLVED)
        ├── v2 (UNRESOLVED)
        │   └── v21 (MATCH)
        │       ├── v211 (UNRESOLVED)
        │       └── v212 (UNRESOLVED)
        └── v3 (UNRESOLVED)
        
        Expected after propagation:
        root (CONSISTENT)
        ├── v1 (CONSISTENT)
        │   ├── v11 (MATCH)         # Preserved
        │   │   └── v111 (MATCH)    # Inherits from MATCH parent
        │   └── v12 (CONSISTENT)    # Inherits from CONSISTENT parent
        ├── v2 (CONSISTENT)         # Inherits from CONSISTENT root
        │   └── v21 (MATCH)         # Preserved
        │       ├── v211 (MATCH)    # Inherits from MATCH parent
        │       └── v212 (MATCH)    # Inherits from MATCH parent
        └── v3 (CONSISTENT)         # Inherits from CONSISTENT root
        """
        # Arrange
        statuses = pd.Series({
            "root": CONSISTENT,
            "v1": CONSISTENT,
            "v11": MATCH,
            "v111": UNRESOLVED,
            "v12": UNRESOLVED,
            "v2": UNRESOLVED,
            "v21": MATCH,
            "v211": UNRESOLVED,
            "v212": UNRESOLVED,
            "v3": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(complex_tree, statuses.copy())
        
        # Assert
        assert result["v11"] == MATCH  # Preserved
        assert result["v111"] == MATCH  # Inherits from MATCH parent
        assert result["v12"] == CONSISTENT  # Inherits from CONSISTENT parent
        assert result["v21"] == MATCH  # Preserved
        assert result["v211"] == MATCH  # Inherits from MATCH parent
        assert result["v212"] == MATCH  # Inherits from MATCH parent
        assert result["v3"] == CONSISTENT  # Inherits from CONSISTENT root

    def test_unclear_propagation(self, simple_tree: dict[str, set[str]]) -> None:
        """Test that UNCLEAR status only affects UNRESOLVED and MAYBE_INCONSISTENT.
        
        Initial state:
        root (UNCLEAR)
        ├── v1 (MATCH)
        │   ├── v11 (UNRESOLVED)
        │   └── v12 (UNRESOLVED)
        └── v2 (INCONSISTENT)
            └── v21 (MAYBE_INCONSISTENT)
        
        Expected after propagation:
        root (UNCLEAR)
        ├── v1 (MATCH)              # Unchanged
        │   ├── v11 (MATCH)         # Inherits from MATCH parent (not UNCLEAR root)
        │   └── v12 (MATCH)         # Inherits from MATCH parent (not UNCLEAR root)
        └── v2 (INCONSISTENT)       # Unchanged
            └── v21 (INCONSISTENT)  # Inherits from INCONSISTENT parent first
        """
        # Arrange
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": MATCH,
            "v11": UNRESOLVED,
            "v12": UNRESOLVED,
            "v2": INCONSISTENT,
            "v21": MAYBE_INCONSISTENT,
        })
        
        # Act
        result = mark_descendant_infections(simple_tree, statuses.copy())
        
        # Assert
        assert result["root"] == UNCLEAR
        assert result["v1"] == MATCH
        assert result["v11"] == MATCH  # Inherits from MATCH parent, not UNCLEAR root
        assert result["v12"] == MATCH  # Inherits from MATCH parent, not UNCLEAR root
        assert result["v2"] == INCONSISTENT
        assert result["v21"] == INCONSISTENT  # Inherits from INCONSISTENT parent first

    def test_unclear_affects_unresolved_and_maybe_inconsistent(self, complex_tree: dict[str, set[str]]) -> None:
        """Test UNCLEAR propagation to UNRESOLVED and MAYBE_INCONSISTENT nodes.
        
        Initial state:
        root (UNCLEAR)
        ├── v1 (UNRESOLVED)
        │   ├── v11 (UNRESOLVED)
        │   │   └── v111 (UNRESOLVED)
        │   └── v12 (MAYBE_INCONSISTENT)
        ├── v2 (MATCH)
        │   └── v21 (CONSISTENT)
        │       ├── v211 (UNRESOLVED)
        │       └── v212 (UNRESOLVED)
        └── v3 (INCONSISTENT)
        
        Expected after propagation:
        root (UNCLEAR)
        ├── v1 (UNCLEAR)            # Changed from UNRESOLVED
        │   ├── v11 (UNCLEAR)       # Changed from UNRESOLVED
        │   │   └── v111 (UNCLEAR)  # Changed from UNRESOLVED
        │   └── v12 (UNCLEAR)       # Changed from MAYBE_INCONSISTENT
        ├── v2 (MATCH)              # Unchanged
        │   └── v21 (MATCH)         # Inherits from MATCH parent (priority)
        │       ├── v211 (MATCH)    # Inherits from MATCH grandparent
        │       └── v212 (MATCH)    # Inherits from MATCH grandparent
        └── v3 (INCONSISTENT)       # Unchanged
        """
        # Arrange
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": UNRESOLVED,
            "v11": UNRESOLVED,
            "v111": UNRESOLVED,
            "v12": MAYBE_INCONSISTENT,
            "v2": MATCH,
            "v21": CONSISTENT,
            "v211": UNRESOLVED,
            "v212": UNRESOLVED,
            "v3": INCONSISTENT,
        })
        
        # Act
        result = mark_descendant_infections(complex_tree, statuses.copy())
        
        # Assert
        assert result["root"] == UNCLEAR
        assert result["v1"] == UNCLEAR  # Changed from UNRESOLVED
        assert result["v11"] == UNCLEAR  # Changed from UNRESOLVED
        assert result["v111"] == UNCLEAR  # Changed from UNRESOLVED
        assert result["v12"] == UNCLEAR  # Changed from MAYBE_INCONSISTENT
        assert result["v2"] == MATCH  # Unchanged
        assert result["v21"] == MATCH  # Inherits from MATCH parent (MATCH propagation happens first)
        assert result["v211"] == MATCH  # Inherits from MATCH parent
        assert result["v212"] == MATCH  # Inherits from MATCH parent
        assert result["v3"] == INCONSISTENT  # Unchanged

    def test_propagation_order_matters(self, simple_tree: dict[str, set[str]]) -> None:
        """Test that propagation order (MATCH > INCONSISTENT > CONSISTENT > UNCLEAR) is respected.
        
        Initial state:
        root (UNCLEAR)
        ├── v1 (MATCH)
        │   ├── v11 (INCONSISTENT)   # Will be overridden
        │   └── v12 (UNRESOLVED)
        └── v2 (INCONSISTENT)
            └── v21 (UNRESOLVED)
        
        Expected after propagation:
        root (UNCLEAR)
        ├── v1 (MATCH)
        │   ├── v11 (MATCH)          # MATCH overrides initial INCONSISTENT
        │   └── v12 (MATCH)          # Inherits from MATCH parent
        └── v2 (INCONSISTENT)
            └── v21 (INCONSISTENT)   # Inherits from INCONSISTENT parent
        """
        # Arrange - Create a situation where order matters
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": MATCH,
            "v11": INCONSISTENT,  # Should become MATCH due to parent
            "v12": UNRESOLVED,
            "v2": INCONSISTENT,
            "v21": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(simple_tree, statuses.copy())
        
        # Assert
        assert result["v1"] == MATCH
        assert result["v11"] == MATCH  # MATCH overrides initial INCONSISTENT
        assert result["v12"] == MATCH
        assert result["v2"] == INCONSISTENT
        assert result["v21"] == INCONSISTENT

    def test_no_unresolved_after_processing(self, complex_tree: dict[str, set[str]]) -> None:
        """Test that no UNRESOLVED statuses remain after processing."""
        # Arrange - Start with many UNRESOLVED nodes
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": UNRESOLVED,
            "v11": UNRESOLVED,
            "v111": UNRESOLVED,
            "v12": UNRESOLVED,
            "v2": UNRESOLVED,
            "v21": UNRESOLVED,
            "v211": UNRESOLVED,
            "v212": UNRESOLVED,
            "v3": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(complex_tree, statuses.copy())
        
        # Assert - All UNRESOLVED should be converted to UNCLEAR
        assert UNRESOLVED not in result.values
        assert all(result == UNCLEAR)

    def test_no_maybe_inconsistent_after_processing(self, complex_tree: dict[str, set[str]]) -> None:
        """Test that no MAYBE_INCONSISTENT statuses remain after processing."""
        # Arrange
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": MAYBE_INCONSISTENT,
            "v11": MAYBE_INCONSISTENT,
            "v111": MAYBE_INCONSISTENT,
            "v12": INCONSISTENT,
            "v2": CONSISTENT,
            "v21": MATCH,
            "v211": MAYBE_INCONSISTENT,
            "v212": MAYBE_INCONSISTENT,
            "v3": MAYBE_INCONSISTENT,
        })
        
        # Act
        result = mark_descendant_infections(complex_tree, statuses.copy())
        
        # Assert - No MAYBE_INCONSISTENT should remain
        assert MAYBE_INCONSISTENT not in result.values

    def test_assertion_unresolved_remaining(self, simple_tree: dict[str, set[str]]) -> None:
        """Test that assertion fires if UNRESOLVED remains (shouldn't happen in practice)."""
        # Arrange - Create a tree with an isolated UNRESOLVED node
        isolated_tree = {"isolated": set()}
        statuses = pd.Series({"isolated": UNRESOLVED})
        
        # Act & Assert
        with pytest.raises(AssertionError, match="Some taxids are still unresolved"):
            mark_descendant_infections(isolated_tree, statuses)

    def test_assertion_maybe_inconsistent_remaining(self, simple_tree: dict[str, set[str]]) -> None:
        """Test that assertion fires if MAYBE_INCONSISTENT remains (shouldn't happen in practice)."""
        # Arrange - Create a tree with an isolated MAYBE_INCONSISTENT node
        isolated_tree = {"isolated": set()}
        statuses = pd.Series({"isolated": MAYBE_INCONSISTENT})
        
        # Act & Assert
        with pytest.raises(AssertionError, match="Some taxids are still maybe inconsistent"):
            mark_descendant_infections(isolated_tree, statuses)

    def test_assertion_consistent_with_inconsistent_descendant(self) -> None:
        """Test that assertion fires if CONSISTENT has INCONSISTENT descendants."""
        # Arrange - This shouldn't happen after mark_ancestor_infections
        tree = {
            "root": {"child"},
            "child": set()
        }
        statuses = pd.Series({
            "root": CONSISTENT,
            "child": INCONSISTENT  # Invalid: CONSISTENT parent can't have INCONSISTENT child
        })
        
        # Act & Assert
        with pytest.raises(AssertionError, match="Some taxids marked as CONSISTENT have descendants marked as INCONSISTENT or UNCLEAR"):
            mark_descendant_infections(tree, statuses)

    def test_assertion_consistent_with_unclear_descendant(self) -> None:
        """Test that assertion fires if CONSISTENT has UNCLEAR descendants."""
        # Arrange
        tree = {
            "root": {"child"},
            "child": set()
        }
        statuses = pd.Series({
            "root": CONSISTENT,
            "child": UNCLEAR  # Invalid: CONSISTENT parent can't have UNCLEAR child
        })
        
        # Act & Assert
        with pytest.raises(AssertionError, match="Some taxids marked as CONSISTENT have descendants marked as INCONSISTENT or UNCLEAR"):
            mark_descendant_infections(tree, statuses)

    def test_realistic_example_from_docs(self, complex_tree: dict[str, set[str]]) -> None:
        """Test a realistic example from cases.md documentation.
        
        Initial state (Example 3 from cases.md):
        root (UNCLEAR)
        ├── v1 (MATCH)
        │   └── v11 (UNRESOLVED)
        ├── v2 (INCONSISTENT)
        │   └── v21 (UNRESOLVED)
        └── v3 (MAYBE_INCONSISTENT)
        
        Expected after propagation:
        root (UNCLEAR)
        ├── v1 (MATCH)
        │   └── v11 (MATCH)         # Inherits from MATCH parent
        ├── v2 (INCONSISTENT)
        │   └── v21 (INCONSISTENT)  # Inherits from INCONSISTENT parent
        └── v3 (UNCLEAR)            # Changed from MAYBE_INCONSISTENT by UNCLEAR parent
        """
        # Arrange - Example 3 from cases.md
        tree = {
            "root": {"v1", "v2", "v3"},
            "v1": {"v11"},
            "v2": {"v21"},
            "v3": set(),
            "v11": set(),
            "v21": set(),
        }
        
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": MATCH,
            "v11": UNRESOLVED,
            "v2": INCONSISTENT,
            "v21": UNRESOLVED,
            "v3": MAYBE_INCONSISTENT,
        })
        
        # Act
        result = mark_descendant_infections(tree, statuses.copy())
        
        # Assert - Should match expected output from documentation
        assert result["root"] == UNCLEAR
        assert result["v1"] == MATCH
        assert result["v11"] == MATCH  # Inherits from MATCH
        assert result["v2"] == INCONSISTENT
        assert result["v21"] == INCONSISTENT  # Inherits from INCONSISTENT
        assert result["v3"] == UNCLEAR  # Changed from MAYBE_INCONSISTENT by UNCLEAR parent

    def test_deep_hierarchy_propagation(self) -> None:
        """Test propagation through a deep hierarchy."""
        # Arrange - Create a deep tree
        tree = {
            "root": {"level1"},
            "level1": {"level2"},
            "level2": {"level3"},
            "level3": {"level4"},
            "level4": {"level5"},
            "level5": set(),
        }
        
        statuses = pd.Series({
            "root": CONSISTENT,
            "level1": UNRESOLVED,
            "level2": UNRESOLVED,
            "level3": UNRESOLVED,
            "level4": UNRESOLVED,
            "level5": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(tree, statuses.copy())
        
        # Assert - All should inherit CONSISTENT
        assert all(result == CONSISTENT)

    def test_complex_mixed_propagation(self) -> None:
        """Test complex scenario with all propagation rules in effect.
        
        Initial state:
        root (UNCLEAR)
        ├── v1 (MATCH)
        │   ├── v11 (INCONSISTENT)      # Will be overridden by MATCH
        │   └── v12 (UNRESOLVED)
        ├── v2 (INCONSISTENT)
        │   ├── v21 (MATCH)             # Will be overridden by INCONSISTENT
        │   └── v22 (UNRESOLVED)
        ├── v3 (CONSISTENT)
        │   └── v31 (UNRESOLVED)
        │       └── v311 (MATCH)        # Will preserve MATCH
        └── v4 (UNRESOLVED)
            ├── v41 (MAYBE_INCONSISTENT)
            └── v42 (UNRESOLVED)
        
        Expected after propagation shows the priority order in action:
        - MATCH propagation (highest priority)
        - INCONSISTENT propagation (second priority)
        - CONSISTENT propagation (third priority, preserves MATCH)
        - UNCLEAR propagation (lowest priority, only affects UNRESOLVED/MAYBE_INCONSISTENT)
        """
        # Arrange
        tree = {
            "root": {"v1", "v2", "v3", "v4"},
            "v1": {"v11", "v12"},
            "v2": {"v21", "v22"},
            "v3": {"v31"},
            "v4": {"v41", "v42"},
            "v11": set(),
            "v12": set(),
            "v21": set(),
            "v22": set(),
            "v31": {"v311"},
            "v311": set(),
            "v41": set(),
            "v42": set(),
        }
        
        statuses = pd.Series({
            "root": UNCLEAR,
            "v1": MATCH,
            "v11": INCONSISTENT,
            "v12": UNRESOLVED,
            "v2": INCONSISTENT,
            "v21": MATCH,  # Should become INCONSISTENT
            "v22": UNRESOLVED,
            "v3": CONSISTENT,
            "v31": UNRESOLVED,
            "v311": MATCH,  # Should stay MATCH
            "v4": UNRESOLVED,
            "v41": MAYBE_INCONSISTENT,
            "v42": UNRESOLVED,
        })
        
        # Act
        result = mark_descendant_infections(tree, statuses.copy())
        
        # Assert
        assert result["v1"] == MATCH
        assert result["v11"] == MATCH  # Overridden by MATCH parent
        assert result["v12"] == MATCH  # Inherits from MATCH parent
        assert result["v2"] == INCONSISTENT
        assert result["v21"] == INCONSISTENT  # Overridden by INCONSISTENT parent
        assert result["v22"] == INCONSISTENT  # Inherits from INCONSISTENT parent
        assert result["v3"] == CONSISTENT
        assert result["v31"] == CONSISTENT  # Inherits from CONSISTENT parent
        assert result["v311"] == MATCH  # Preserved (CONSISTENT doesn't override MATCH)
        assert result["v4"] == UNCLEAR  # Inherits from UNCLEAR root
        assert result["v41"] == UNCLEAR  # Changed from MAYBE_INCONSISTENT
        assert result["v42"] == UNCLEAR  # Changed from UNRESOLVED


# =======================================================================
# Main execution
# =======================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
