"""Unit tests for domain segmentation and regional identity identification."""

import sqlite3

import pandas as pd
import pytest

import numpy as np
from src.beglib import (
    calculate_shannon_entropy,
    collapse_to_canonical,
    identify_regional_domains,
    prepare_kmer_frequencies,
    write_domains_db,
)


@pytest.fixture
def mock_kmer_df() -> pd.DataFrame:
    """Fixture providing synthetic genomic window data with 2 distinct k-mer regions and an N-gap."""
    data = {
        "library_name": ["lib1"] * 10,
        "seq_name": ["chr1"] * 10,
        "window_start": list(range(1, 101, 10)),
        "window_stop": list(range(10, 101, 10)),
        "window_size": [10] * 10,
        # First 4 windows AT-rich
        "AAA": [5, 4, 5, 4, 0, 1, 0, 1, 0, 1],
        "TTT": [4, 5, 4, 5, 0, 1, 0, 1, 0, 1],
        # 5th window is an N-gap (all counts 0)
        "CCC": [0, 0, 0, 0, 0, 4, 5, 4, 5, 4],
        "GGG": [0, 0, 0, 0, 0, 4, 4, 5, 4, 5],
        # Remaining columns
        "AAC": [1, 1, 1, 1, 0, 0, 0, 0, 1, 0],
        "GTT": [1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    }
    return pd.DataFrame(data)


def test_collapse_to_canonical(mock_kmer_df: pd.DataFrame) -> None:
    """Test collapse_to_canonical combines reverse complement k-mers."""
    kmer_cols = ["AAA", "TTT", "CCC", "GGG", "AAC", "GTT"]
    canon_df = collapse_to_canonical(mock_kmer_df, kmer_cols)

    assert "TTT" not in canon_df.columns
    assert "AAA" in canon_df.columns
    # AAA + TTT for row 0: 5 + 4 = 9
    assert canon_df["AAA"].iloc[0] == 9
    # CCC + GGG for row 5: 4 + 4 = 8
    assert canon_df["CCC"].iloc[5] == 8


def test_calculate_shannon_entropy() -> None:
    """Test calculate_shannon_entropy produces expected natural entropy values."""
    df_equal = pd.DataFrame({"k1": [0.25], "k2": [0.25], "k3": [0.25], "k4": [0.25]})
    entropy_equal = calculate_shannon_entropy(df_equal)
    assert pytest.approx(entropy_equal.iloc[0], abs=1e-4) == np.log(4)

    df_pure = pd.DataFrame({"k1": [1.0], "k2": [0.0], "k3": [0.0]})
    entropy_pure = calculate_shannon_entropy(df_pure)
    assert pytest.approx(entropy_pure.iloc[0]) == 0.0

    df_zero = pd.DataFrame({"k1": [0.0], "k2": [0.0]})
    entropy_zero = calculate_shannon_entropy(df_zero)
    assert pytest.approx(entropy_zero.iloc[0]) == 0.0


def test_identify_regional_domains_with_smoothing_and_complexity(mock_kmer_df: pd.DataFrame) -> None:
    """Test identify_regional_domains with custom smoothing sigma and complexity column."""
    mock_df_with_comp = mock_kmer_df.copy()
    mock_df_with_comp["complexity"] = [0.5, 0.52, 0.51, 0.49, 0.0, 0.8, 0.82, 0.81, 0.79, 0.8]

    kmer_cols = ["AAA", "TTT", "CCC", "GGG", "AAC", "GTT"]
    df_results, domain_profiles = identify_regional_domains(
        df=mock_df_with_comp,
        kmer_cols=kmer_cols,
        n_states="2",
        smoothing_sigma=1.5,
    )

    assert "Expected_2" in df_results.columns
    assert "complexity" in df_results.columns
    assert len(df_results) == 10
    assert len(domain_profiles) == 2


def test_prepare_kmer_frequencies(mock_kmer_df: pd.DataFrame) -> None:
    """Test prepare_kmer_frequencies normalizes counts to proportions."""
    kmer_cols = ["AAA", "TTT", "CCC", "GGG", "AAC", "GTT"]
    freq_df, features = prepare_kmer_frequencies(mock_kmer_df, kmer_cols, use_canonical=False)

    assert len(features) == len(kmer_cols)
    assert freq_df.iloc[0].sum() == pytest.approx(1.0)
    # Zero row should yield all zeros without division by zero errors
    assert freq_df.iloc[4].sum() == pytest.approx(0.0)


def test_identify_regional_domains(mock_kmer_df: pd.DataFrame) -> None:
    """Test identify_regional_domains segments windows into domain states."""
    kmer_cols = ["AAA", "TTT", "CCC", "GGG", "AAC", "GTT"]
    df_results, domain_profiles = identify_regional_domains(
        df=mock_kmer_df,
        kmer_cols=kmer_cols,
        use_canonical=False,
        n_states="2",
        variance_target=0.90,
    )

    assert "Expected_2" in df_results.columns
    assert "comment" in df_results.columns
    assert "entropy" in df_results.columns
    assert "complexity" in df_results.columns
    assert "AAA" not in df_results.columns  # Truncated output check

    # Verify N-gap comment for 5th window (index 4)
    assert df_results["comment"].iloc[4] == "N-gap present"
    assert df_results["comment"].iloc[0] == ""

    # Verify domain profiles summary table structure
    assert "Expected_States" in domain_profiles.columns
    assert "Domain_ID" in domain_profiles.columns
    assert "entropy" in domain_profiles.columns
    assert "complexity" in domain_profiles.columns
    assert len(domain_profiles) == 2
    assert set(kmer_cols).issubset(set(domain_profiles.columns))


def test_identify_regional_domains_range(mock_kmer_df: pd.DataFrame) -> None:
    """Test identify_regional_domains with range of states (e.g. '2:3')."""
    kmer_cols = ["AAA", "TTT", "CCC", "GGG", "AAC", "GTT"]
    df_results, domain_profiles = identify_regional_domains(
        df=mock_kmer_df,
        kmer_cols=kmer_cols,
        n_states="1:2",
    )

    assert "Expected_1" in df_results.columns
    assert "Expected_2" in df_results.columns
    assert len(domain_profiles) == 3  # 1 for state count 1, 2 for state count 2


def test_write_domains_db(tmp_path, mock_kmer_df: pd.DataFrame) -> None:
    """Test writing domain predictions and profiles into SQLite tables."""
    kmer_cols = ["AAA", "TTT", "CCC", "GGG", "AAC", "GTT"]
    df_results, domain_profiles = identify_regional_domains(
        df=mock_kmer_df,
        kmer_cols=kmer_cols,
        n_states="2",
    )

    db_path = tmp_path / "domains.db"
    conn = sqlite3.connect(str(db_path))

    write_domains_db(
        conn=conn,
        df_results=df_results,
        domain_profiles=domain_profiles,
        params={"n_states": "2", "variance_target": 0.90, "canonical": False},
    )

    cursor = conn.cursor()
    cursor.execute("SELECT run_id, expected_states FROM domain_runs")
    runs = cursor.fetchall()
    assert len(runs) == 1
    run_id = runs[0][0]
    assert run_id == "run_n2_v0.9_c0"

    cursor.execute("SELECT run_id, seq_name, domain_id, comment FROM domains WHERE run_id = ?", (run_id,))
    domain_rows = cursor.fetchall()
    assert len(domain_rows) == 10
    assert domain_rows[4][3] == "N-gap present"

    cursor.execute("SELECT run_id, domain_id, profile_json FROM domain_profiles WHERE run_id = ?", (run_id,))
    profile_rows = cursor.fetchall()
    assert len(profile_rows) == 2

    conn.close()
