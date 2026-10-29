"""Tests for neighbor block collation in domain identification."""

import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import sqlite3

import src.beglib as bl
from beg import build_parser
from src.utils.domain_visualize import visualize_domains


def test_collate_neighboring_blocks_identical():
    """Test that identical neighboring windows are collated into a single block."""
    df_results = pd.DataFrame({
        "library_name": ["lib1", "lib1", "lib1"],
        "seq_name": ["chr1", "chr1", "chr1"],
        "window_start": [1, 10001, 20001],
        "window_stop": [10000, 20000, 30000],
        "window_size": [10000, 10000, 10000],
        "complexity": [0.5, 0.5, 0.5],
        "entropy": [2.0, 2.0, 2.0],
        "comment": ["", "", ""],
    })

    kmer_cols = ["AAA", "AAC", "AAG"]
    model_freq_df = pd.DataFrame({
        "AAA": [0.5, 0.5, 0.5],
        "AAC": [0.3, 0.3, 0.3],
        "AAG": [0.2, 0.2, 0.2],
        "entropy": [2.0, 2.0, 2.0],
        "complexity": [0.5, 0.5, 0.5],
    })

    collated_res, collated_freq = bl.collate_neighboring_blocks(
        df_results=df_results,
        model_freq_df=model_freq_df,
        feature_cols=kmer_cols,
        collation_threshold=0.005,
    )

    assert len(collated_res) == 1
    assert collated_res["window_start"].iloc[0] == 1
    assert collated_res["window_stop"].iloc[0] == 30000
    assert collated_res["window_size"].iloc[0] == 30000
    assert pytest.approx(collated_res["complexity"].iloc[0]) == 0.5
    assert pytest.approx(collated_res["entropy"].iloc[0]) == 2.0
    assert pytest.approx(collated_freq["AAA"].iloc[0]) == 0.5


def test_collate_neighboring_blocks_different_contigs():
    """Test that windows from different sequence contigs are never collated."""
    df_results = pd.DataFrame({
        "library_name": ["lib1", "lib1"],
        "seq_name": ["chr1", "chr2"],
        "window_start": [1, 1],
        "window_stop": [10000, 10000],
        "window_size": [10000, 10000],
        "complexity": [0.5, 0.5],
        "entropy": [2.0, 2.0],
        "comment": ["", ""],
    })

    kmer_cols = ["AAA", "AAC"]
    model_freq_df = pd.DataFrame({
        "AAA": [0.5, 0.5],
        "AAC": [0.5, 0.5],
        "entropy": [2.0, 2.0],
        "complexity": [0.5, 0.5],
    })

    collated_res, collated_freq = bl.collate_neighboring_blocks(
        df_results=df_results,
        model_freq_df=model_freq_df,
        feature_cols=kmer_cols,
        collation_threshold=0.05,
    )

    assert len(collated_res) == 2
    assert list(collated_res["seq_name"]) == ["chr1", "chr2"]


def test_collate_neighboring_blocks_weighted_averaging():
    """Test length-weighted averaging of features when merging unequal window sizes."""
    df_results = pd.DataFrame({
        "library_name": ["lib1", "lib1"],
        "seq_name": ["chr1", "chr1"],
        "window_start": [1, 10001],
        "window_stop": [10000, 30000],
        "window_size": [10000, 20000],  # Weight ratio 1 : 2
        "complexity": [0.3, 0.6],
        "entropy": [1.0, 2.0],
        "comment": ["", ""],
    })

    kmer_cols = ["AAA"]
    model_freq_df = pd.DataFrame({
        "AAA": [0.3, 0.6],
        "entropy": [1.0, 2.0],
        "complexity": [0.3, 0.6],
    })

    # High threshold to force merge
    collated_res, collated_freq = bl.collate_neighboring_blocks(
        df_results=df_results,
        model_freq_df=model_freq_df,
        feature_cols=kmer_cols,
        collation_threshold=1.0,
    )

    assert len(collated_res) == 1
    # Expected weighted average: (0.3*10000 + 0.6*20000) / 30000 = 0.5
    assert pytest.approx(collated_res["complexity"].iloc[0]) == 0.5
    expected_entropy = (1.0 * 10000 + 2.0 * 20000) / 30000
    assert pytest.approx(collated_res["entropy"].iloc[0]) == expected_entropy
    assert pytest.approx(collated_freq["entropy"].iloc[0]) == expected_entropy


def test_cli_parser_collation_threshold():
    """Test that CLI parser correctly parses --collation-threshold."""
    parser = build_parser()
    args = parser.parse_args(["domains", "-i", "input.beg", "--collation-threshold", "0.02"])
    assert args.collation_threshold == 0.02


def test_identify_regional_domains_with_collation():
    """Test end-to-end domain identification with collation threshold."""
    np.random.seed(42)
    kmer_cols = [f"kmer_{i}" for i in range(16)]
    rows = []

    # 40 total windows: 20 in block 1, 20 in block 2
    for i in range(40):
        if i < 20:
            counts = np.random.poisson(lam=10, size=16)
        else:
            counts = np.random.poisson(lam=100, size=16)

        row = {
            "library_name": "lib",
            "seq_name": "chr1",
            "window_start": i * 1000 + 1,
            "window_stop": (i + 1) * 1000,
            "window_size": 1000,
            "complexity": 0.5 + (0.01 if i >= 20 else 0.0),
        }
        for idx, k in enumerate(kmer_cols):
            row[k] = counts[idx]
        rows.append(row)

    df = pd.DataFrame(rows)

    df_results, profiles = bl.identify_regional_domains(
        df=df,
        kmer_cols=kmer_cols,
        n_states=2,
        collation_threshold=0.08,
        smoothing_sigma=0.0,
    )

    assert len(df_results) < len(df)
    assert "Expected_2" in df_results.columns
    assert not profiles.empty


def test_visualize_domains_collated():
    """Test visualize_domains utility with collated domain blocks."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        counts_file = tmp_path / "counts.beg"
        domains_file = tmp_path / "domains.beg"

        # 3 raw counts windows
        counts_df = pd.DataFrame({
            "library_name": ["lib", "lib", "lib"],
            "seq_name": ["chr1", "chr1", "chr1"],
            "window_start": [1, 10001, 20001],
            "window_stop": [10000, 20000, 30000],
            "window_size": [10000, 10000, 10000],
            "complexity": [0.5, 0.5, 0.5],
            "AAA": [10, 10, 100],
            "AAC": [20, 20, 200],
        })
        counts_df.to_csv(counts_file, sep="\t", index=False)

        # 2 collated domain blocks
        domains_df = pd.DataFrame({
            "library_name": ["lib", "lib"],
            "seq_name": ["chr1", "chr1"],
            "window_start": [1, 20001],
            "window_stop": [20000, 30000],
            "window_size": [20000, 10000],
            "complexity": [0.5, 0.5],
            "entropy": [1.5, 2.5],
            "Expected_2": ["Domain_0", "Domain_1"],
        })
        domains_df.to_csv(domains_file, sep="\t", index=False)

        visualize_domains(
            domains_path=domains_file,
            counts_path=counts_file,
            output_dir=tmp_path / "plots",
            img_format="png",
        )

        out_img = tmp_path / "plots" / "chr1_domains.png"
        assert out_img.exists()
