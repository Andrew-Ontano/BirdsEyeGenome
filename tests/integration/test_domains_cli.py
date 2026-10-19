"""Integration tests for the `beg domains` CLI subcommand."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from src.beglib import init_db, write_counts_beg, write_counts_db


@pytest.fixture
def sample_beg_file(tmp_path: Path) -> Path:
    """Fixture to generate a sample .beg file for CLI testing."""
    beg_path = tmp_path / "sample.beg"
    kmers = ["AAA", "TTT", "CCC", "GGG"]
    results = [
        {"seq_name": "chr1", "start": 0, "end": 10, "window_size": 10, "counts": {"AAA": 5, "TTT": 5, "CCC": 0, "GGG": 0}},
        {"seq_name": "chr1", "start": 10, "end": 20, "window_size": 10, "counts": {"AAA": 4, "TTT": 6, "CCC": 0, "GGG": 0}},
        {"seq_name": "chr1", "start": 20, "end": 30, "window_size": 10, "counts": {"AAA": 0, "TTT": 0, "CCC": 5, "GGG": 5}},
        {"seq_name": "chr1", "start": 30, "end": 40, "window_size": 10, "counts": {"AAA": 0, "TTT": 0, "CCC": 6, "GGG": 4}},
    ]
    write_counts_beg(str(beg_path), "sample.fa", "sample.fa", kmers, results)
    return beg_path


@pytest.fixture
def sample_db_file(tmp_path: Path) -> Path:
    """Fixture to generate a sample .db file for CLI testing."""
    db_path = tmp_path / "sample.db"
    kmers = ["AAA", "TTT", "CCC", "GGG"]
    conn = init_db(str(db_path), kmers)
    results = [
        {"seq_name": "chr1", "start": 0, "end": 10, "window_size": 10, "counts": {"AAA": 5, "TTT": 5, "CCC": 0, "GGG": 0}},
        {"seq_name": "chr1", "start": 10, "end": 20, "window_size": 10, "counts": {"AAA": 4, "TTT": 6, "CCC": 0, "GGG": 0}},
        {"seq_name": "chr1", "start": 20, "end": 30, "window_size": 10, "counts": {"AAA": 0, "TTT": 0, "CCC": 5, "GGG": 5}},
        {"seq_name": "chr1", "start": 30, "end": 40, "window_size": 10, "counts": {"AAA": 0, "TTT": 0, "CCC": 6, "GGG": 4}},
    ]
    write_counts_db(conn, "sample.fa", results, kmers)
    conn.close()
    return db_path


def test_beg_domains_cli_beg_input(sample_beg_file: Path, tmp_path: Path) -> None:
    """Test `beg domains` CLI with a .beg input file."""
    out_beg = tmp_path / "out_domains.beg"
    out_summary = tmp_path / "out_profiles.tsv"

    cmd = [
        sys.executable,
        "beg.py",
        "domains",
        "-i",
        str(sample_beg_file),
        "-o",
        str(out_beg),
        "-s",
        str(out_summary),
        "-n",
        "2",
        "--canonical",
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"CLI error: {res.stderr}"

    assert out_beg.exists()
    assert out_summary.exists()

    df_out = pd.read_csv(out_beg, sep="\t")
    assert "Expected_2" in df_out.columns
    assert "comment" in df_out.columns
    assert "entropy" in df_out.columns
    assert "complexity" in df_out.columns
    assert len(df_out) == 4

    df_sum = pd.read_csv(out_summary, sep="\t")
    assert "Expected_States" in df_sum.columns
    assert "Domain_ID" in df_sum.columns
    assert "entropy" in df_sum.columns
    assert "complexity" in df_sum.columns
    assert len(df_sum) == 2


def test_beg_domains_cli_smoothing_sigma(sample_beg_file: Path, tmp_path: Path) -> None:
    """Test `beg domains` CLI with --smoothing-sigma option."""
    out_beg = tmp_path / "out_smoothed_domains.beg"
    cmd = [
        sys.executable,
        "beg.py",
        "domains",
        "-i",
        str(sample_beg_file),
        "-o",
        str(out_beg),
        "-n",
        "2",
        "--smoothing-sigma",
        "1.5",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"CLI error: {res.stderr}"
    assert out_beg.exists()
    df_out = pd.read_csv(out_beg, sep="\t")
    assert "Expected_2" in df_out.columns
    assert "complexity" in df_out.columns


def test_beg_domains_cli_db_input(sample_db_file: Path, tmp_path: Path) -> None:
    """Test `beg domains` CLI with a .db input file."""
    out_db = tmp_path / "out_domains.db"
    out_summary = tmp_path / "out_profiles.tsv"

    cmd = [
        sys.executable,
        "beg.py",
        "domains",
        "-i",
        str(sample_db_file),
        "-o",
        str(out_db),
        "-s",
        str(out_summary),
        "-n",
        "2",
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"CLI error: {res.stderr}"

    assert out_db.exists()
    assert out_summary.exists()
