"""Unit tests for src/utils/hierarchical_composition.py."""

import tempfile
from pathlib import Path
import pandas as pd
import pytest
import sqlite3

from src.utils.hierarchical_composition import (
    load_domain_dataset,
    map_parent_domain,
    build_hierarchical_composition,
)


@pytest.fixture
def sample_domain_files(tmp_path):
    """Create temporary TSV domain files representing Small (1000bp) and Large (2000bp) window sizes."""
    # Small window dataset (1000 bp)
    small_df = pd.DataFrame({
        "library_name": ["lib1", "lib1"],
        "seq_name": ["chr1", "chr1"],
        "window_start": [1, 1001],
        "window_stop": [1000, 2000],
        "window_size": [1000, 1000],
        "complexity": [0.5, 0.6],
        "entropy": [2.1, 2.3],
        "Expected_2": ["Domain_0", "Domain_1"],
        "Expected_3": ["Domain_0", "Domain_2"],
    })
    small_path = tmp_path / "small_domains.beg"
    small_df.to_csv(small_path, sep="\t", index=False)

    # Large window dataset (2000 bp)
    large_df = pd.DataFrame({
        "library_name": ["lib1"],
        "seq_name": ["chr1"],
        "window_start": [1],
        "window_stop": [2000],
        "window_size": [2000],
        "complexity": [0.55],
        "entropy": [2.2],
        "Expected_2": ["Domain_1"],
        "Expected_3": ["Domain_1"],
    })
    large_path = tmp_path / "large_domains.beg"
    large_df.to_csv(large_path, sep="\t", index=False)

    return small_path, large_path


def test_load_domain_dataset(sample_domain_files):
    small_path, large_path = sample_domain_files

    df_s, w_size_s = load_domain_dataset(small_path)
    assert w_size_s == 1000
    assert len(df_s) == 2
    assert "Expected_2" in df_s.columns

    df_l, w_size_l = load_domain_dataset(large_path)
    assert w_size_l == 2000
    assert len(df_l) == 1


def test_map_parent_domain():
    base_df = pd.DataFrame({
        "seq_name": ["chr1", "chr1"],
        "window_start": [1, 1001],
        "window_stop": [1000, 2000],
    })
    target_df = pd.DataFrame({
        "seq_name": ["chr1"],
        "window_start": [1],
        "window_stop": [2000],
        "Expected_2": ["Domain_A"],
    })

    mapped = map_parent_domain(base_df, target_df, "Expected_2")
    assert mapped == ["Domain_A", "Domain_A"]


def test_build_hierarchical_composition(sample_domain_files):
    small_path, large_path = sample_domain_files

    # Pass in unsorted file order
    res_df = build_hierarchical_composition([small_path, large_path])

    # Should have 2 base windows (from small_domains.beg)
    assert len(res_df) == 2
    assert list(res_df["window_size"]) == [1000, 1000]

    # Verify column ordering
    expected_meta = [
        "library_name",
        "seq_name",
        "window_start",
        "window_stop",
        "window_size",
        "complexity",
        "entropy",
        "hierarchical_label",
    ]
    for col in expected_meta:
        assert col in res_df.columns

    # Verify Level columns exist: Expected_2 (Large, Small) -> Level_1, Level_2; Expected_3 (Large, Small) -> Level_3, Level_4
    assert "Level_1" in res_df.columns
    assert "Level_2" in res_df.columns
    assert "Level_3" in res_df.columns
    assert "Level_4" in res_df.columns

    # Check Level_1 (Large Expected_2) mapping -> both windows are in Large Domain_1
    assert list(res_df["Level_1"]) == ["Domain_1", "Domain_1"]
    # Check Level_2 (Small Expected_2) mapping -> Domain_0, Domain_1
    assert list(res_df["Level_2"]) == ["Domain_0", "Domain_1"]

    # Check hierarchical_label shorthand formatting
    # Level_1: 1, Level_2: 0, Level_3: 1, Level_4: 0 -> "1-0-1-0" for row 0
    # Level_1: 1, Level_2: 1, Level_3: 1, Level_4: 2 -> "1-1-1-2" for row 1
    assert res_df["hierarchical_label"].iloc[0] == "1-0-1-0"
    assert res_df["hierarchical_label"].iloc[1] == "1-1-1-2"
