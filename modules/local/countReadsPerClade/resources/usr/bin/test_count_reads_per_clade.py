from collections import Counter

import pytest

from count_reads_per_clade import (
    build_tree,
    children,
    count_direct_reads_per_taxid,
    get_clade_counts,
    is_duplicate,
    parents,
    roots,
)


def test_is_duplicate():
    # Test case where read is not a duplicate (seq_id matches exemplar)
    read_not_duplicate = {"seq_id": "read123", "prim_align_dup_exemplar": "read123"}
    assert not is_duplicate(read_not_duplicate)

    # Test case where read is a duplicate (seq_id differs from exemplar)
    read_is_duplicate = {"seq_id": "read456", "prim_align_dup_exemplar": "read123"}
    assert is_duplicate(read_is_duplicate)

    # Test KeyError when seq_id is missing
    with pytest.raises(KeyError):
        is_duplicate({"prim_align_dup_exemplar": "read123"})

    # Test KeyError when prim_align_dup_exemplar is missing
    with pytest.raises(KeyError):
        is_duplicate({"seq_id": "read123"})

    # Test KeyError when both fields are missing
    with pytest.raises(KeyError):
        is_duplicate({})


def test_parents():
    # Test with a simple tree: 1->2, 1->3, 2->4
    # Parents should be: {1, 2} (nodes that have children)
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "3", "parent_taxid": "1"},
        {"taxid": "4", "parent_taxid": "2"},
    ]
    tree = build_tree(iter(tax_data))
    result = parents(tree)
    expected = {1, 2}
    assert result == expected

    # Test with empty tree
    empty_tree = build_tree(iter([]))
    result = parents(empty_tree)
    assert result == set()

    # Test with single parent-child relationship
    tax_data = [{"taxid": "20", "parent_taxid": "10"}]
    single_tree = build_tree(iter(tax_data))
    result = parents(single_tree)
    assert result == {10}

    # Test return type is set
    tax_data = [{"taxid": "2", "parent_taxid": "1"}]
    tree = build_tree(iter(tax_data))
    result = parents(tree)
    assert isinstance(result, set)


def test_children():
    # Test with a simple tree: 1->2, 1->3, 2->4
    # Children should be: {2, 3, 4} (nodes that are children of other nodes)
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "3", "parent_taxid": "1"},
        {"taxid": "4", "parent_taxid": "2"},
    ]
    tree = build_tree(iter(tax_data))
    result = children(tree)
    expected = {2, 3, 4}
    assert result == expected

    # Test with empty tree
    empty_tree = build_tree(iter([]))
    result = children(empty_tree)
    assert result == set()

    # Test with single parent-child relationship
    tax_data = [{"taxid": "20", "parent_taxid": "10"}]
    single_tree = build_tree(iter(tax_data))
    result = children(single_tree)
    assert result == {20}

    # Test with duplicate children in different parents (should still return unique set)
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "3", "parent_taxid": "1"},
        {"taxid": "2", "parent_taxid": "4"},
        {"taxid": "5", "parent_taxid": "4"},
    ]
    tree_with_duplicates = build_tree(iter(tax_data))
    result = children(tree_with_duplicates)
    expected = {2, 3, 5}
    assert result == expected

    # Test return type is set
    tax_data = [{"taxid": "2", "parent_taxid": "1"}]
    tree = build_tree(iter(tax_data))
    result = children(tree)
    assert isinstance(result, set)


def test_roots():
    # Test with a simple tree: 1->2, 1->3, 2->4
    # Root should be: {1} (parent that is never a child)
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "3", "parent_taxid": "1"},
        {"taxid": "4", "parent_taxid": "2"},
    ]
    tree = build_tree(iter(tax_data))
    result = roots(tree)
    expected = {1}
    assert result == expected

    # Test with empty tree
    empty_tree = build_tree(iter([]))
    result = roots(empty_tree)
    assert result == set()

    # Test with multiple roots (forest): 1->2, 3->4
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "4", "parent_taxid": "3"},
    ]
    forest = build_tree(iter(tax_data))
    result = roots(forest)
    expected = {1, 3}
    assert result == expected

    # Test with single parent-child relationship
    tax_data = [{"taxid": "20", "parent_taxid": "10"}]
    single_tree = build_tree(iter(tax_data))
    result = roots(single_tree)
    assert result == {10}

    # Test with cycle-like structure: 1->2, 2->3, 4->1
    # Root should be {4} (only parent that's never a child)
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "3", "parent_taxid": "2"},
        {"taxid": "1", "parent_taxid": "4"},
    ]
    tree_with_cycle_parent = build_tree(iter(tax_data))
    result = roots(tree_with_cycle_parent)
    expected = {4}
    assert result == expected

    # Test return type is set
    tax_data = [{"taxid": "2", "parent_taxid": "1"}]
    tree = build_tree(iter(tax_data))
    result = roots(tree)
    assert isinstance(result, set)


def test_build_tree():
    # Test basic tree building
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "3", "parent_taxid": "1"},
        {"taxid": "4", "parent_taxid": "2"},
    ]
    result = build_tree(iter(tax_data))
    expected = {1: {2, 3}, 2: {4}}
    assert result == expected

    # Test with custom field names
    tax_data = [
        {"child_id": "10", "parent_id": "5"},
        {"child_id": "11", "parent_id": "5"},
    ]
    result = build_tree(
        iter(tax_data), child_field="child_id", parent_field="parent_id"
    )
    expected = {5: {10, 11}}
    assert result == expected

    # Test with duplicate relationships (should only include each child once)
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "2", "parent_taxid": "1"},  # duplicate
        {"taxid": "3", "parent_taxid": "1"},
    ]
    result = build_tree(iter(tax_data))
    expected = {1: {2, 3}}
    assert result == expected

    # Test with empty data
    result = build_tree(iter([]))
    assert result == {}

    # Test string to int conversion
    tax_data = [{"taxid": "100", "parent_taxid": "50"}]
    result = build_tree(iter(tax_data))
    # Keys and values should be integers
    assert 50 in result
    assert 100 in result[50]


def test_count_direct_reads_per_taxid():
    # Test basic counting with some duplicates
    read_data = [
        {
            "aligner_taxid_lca": "100",
            "seq_id": "read1",
            "prim_align_dup_exemplar": "read1",
        },  # not duplicate
        {
            "aligner_taxid_lca": "100",
            "seq_id": "read2",
            "prim_align_dup_exemplar": "read1",
        },  # duplicate
        {
            "aligner_taxid_lca": "200",
            "seq_id": "read3",
            "prim_align_dup_exemplar": "read3",
        },  # not duplicate
        {
            "aligner_taxid_lca": "100",
            "seq_id": "read4",
            "prim_align_dup_exemplar": "read4",
        },  # not duplicate
    ]

    total, dedup = count_direct_reads_per_taxid(iter(read_data))

    # Total counts: taxid 100 has 3 reads, taxid 200 has 1 read
    assert total[100] == 3
    assert total[200] == 1

    # Deduplicated counts: taxid 100 has 2 non-duplicate reads, taxid 200 has 1
    assert dedup[100] == 2
    assert dedup[200] == 1

    # Test with custom taxid field
    read_data = [
        {"custom_taxid": "50", "seq_id": "read1", "prim_align_dup_exemplar": "read1"}
    ]
    total, dedup = count_direct_reads_per_taxid(
        iter(read_data), taxid_field="custom_taxid"
    )
    assert total[50] == 1
    assert dedup[50] == 1

    # Test with empty data
    total, dedup = count_direct_reads_per_taxid(iter([]))
    assert len(total) == 0
    assert len(dedup) == 0

    # Test return types are Counters
    total, dedup = count_direct_reads_per_taxid(iter([]))
    assert isinstance(total, Counter)
    assert isinstance(dedup, Counter)


def test_get_clade_counts():
    # Test with simple tree: 1->2, 1->3, 2->4
    # If node counts are: {1:5, 2:10, 3:7, 4:3}
    # Then clade counts should be:
    # - Node 4: 3 (only itself)
    # - Node 3: 7 (only itself)
    # - Node 2: 13 (10 + 3 from child 4)
    # - Node 1: 25 (5 + 10 + 3 + 7 from all descendants)
    tax_data = [
        {"taxid": "2", "parent_taxid": "1"},
        {"taxid": "3", "parent_taxid": "1"},
        {"taxid": "4", "parent_taxid": "2"},
    ]
    tree = build_tree(iter(tax_data))
    node_counts = Counter({1: 5, 2: 10, 3: 7, 4: 3})

    result = get_clade_counts(node_counts, tree)
    expected = Counter({1: 25, 2: 13, 3: 7, 4: 3})
    assert result == expected

    # Test with nodes that have zero counts
    tax_data = [{"taxid": "2", "parent_taxid": "1"}]
    tree = build_tree(iter(tax_data))
    node_counts = Counter({1: 5})  # node 2 has 0 count
    result = get_clade_counts(node_counts, tree)
    expected = Counter({1: 5, 2: 0})
    assert result == expected

    # Test with empty tree
    tree = build_tree(iter([]))
    node_counts = Counter({1: 5})  # node 1 not in tree
    result = get_clade_counts(node_counts, tree)
    expected = Counter({})  # empty result since no nodes in tree
    assert result == expected

    # Test with node_counts containing nodes not in tree (should be ignored)
    tax_data = [{"taxid": "2", "parent_taxid": "1"}]
    tree = build_tree(iter(tax_data))
    node_counts = Counter({1: 3, 2: 7, 999: 100})  # 999 not in tree
    result = get_clade_counts(node_counts, tree)
    expected = Counter({1: 10, 2: 7})  # 999 ignored

    # Test return type is Counter
    tax_data = [{"taxid": "2", "parent_taxid": "1"}]
    tree = build_tree(iter(tax_data))
    node_counts = Counter({1: 1})
    result = get_clade_counts(node_counts, tree)
    assert isinstance(result, Counter)
