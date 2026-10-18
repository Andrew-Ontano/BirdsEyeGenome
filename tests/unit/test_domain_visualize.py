"""Unit tests for domain_visualize utility script."""

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from src.beglib import identify_regional_domains, init_db, write_counts_beg, write_counts_db, write_domains_db
from src.utils.domain_visualize import load_domain_data, visualize_domains


@pytest.fixture
def sample_counts_and_domains(tmp_path: Path):
    """Fixture to generate matching count and domain files in .beg and .db formats."""
    kmers = ["AAA", "TTT", "CCC", "GGG"]
    count_results = [
        {"seq_name": "chr1", "start": 0, "end": 10, "window_size": 10, "counts": {"AAA": 5, "TTT": 5, "CCC": 0, "GGG": 0}},
        {"seq_name": "chr1", "start": 10, "end": 20, "window_size": 10, "counts": {"AAA": 4, "TTT": 6, "CCC": 0, "GGG": 0}},
        {"seq_name": "chr1", "start": 20, "end": 30, "window_size": 10, "counts": {"AAA": 0, "TTT": 0, "CCC": 5, "GGG": 5}},
        {"seq_name": "chr1", "start": 30, "end": 40, "window_size": 10, "counts": {"AAA": 0, "TTT": 0, "CCC": 6, "GGG": 4}},
    ]

    # BEG count
    counts_beg = tmp_path / "counts.beg"
    write_counts_beg(str(counts_beg), "sample.fa", "sample.fa", kmers, count_results)

    # BEG domains
    df_counts = pd.DataFrame({
        "library_name": ["sample.fa"] * 4,
        "seq_name": ["chr1"] * 4,
        "window_start": [1, 11, 21, 31],
        "window_stop": [10, 20, 30, 40],
        "window_size": [10] * 4,
        "AAA": [5, 4, 0, 0],
        "TTT": [5, 6, 0, 0],
        "CCC": [0, 0, 5, 6],
        "GGG": [0, 0, 5, 4],
    })
    df_results, domain_profiles = identify_regional_domains(df_counts, kmer_cols=kmers, n_states="2:3")

    domains_beg = tmp_path / "domains.beg"
    df_results.to_csv(domains_beg, sep="\t", index=False)

    # DB counts and domains
    db_file = tmp_path / "sample.db"
    conn = init_db(str(db_file), kmers)
    write_counts_db(conn, "sample.fa", count_results, kmers)
    write_domains_db(conn, df_results, domain_profiles)
    conn.close()

    return counts_beg, domains_beg, db_file


def test_load_domain_data(sample_counts_and_domains) -> None:
    """Test load_domain_data correctly loads expected state columns from .beg and .db files."""
    _counts_beg, domains_beg, db_file = sample_counts_and_domains

    df_beg = load_domain_data(domains_beg)
    assert "Expected_2" in df_beg.columns
    assert "Expected_3" in df_beg.columns

    df_db = load_domain_data(db_file)
    assert "Expected_2" in df_db.columns
    assert "Expected_3" in df_db.columns


def test_visualize_domains_beg(sample_counts_and_domains, tmp_path: Path) -> None:
    """Test visualize_domains creates visualization output image from .beg files."""
    counts_beg, domains_beg, _db_file = sample_counts_and_domains
    out_dir = tmp_path / "viz_beg_out"

    visualize_domains(
        domains_path=domains_beg,
        counts_path=counts_beg,
        output_dir=out_dir,
    )

    expected_img = out_dir / "chr1_domains.png"
    assert expected_img.exists()
    assert expected_img.stat().st_size > 0


def test_visualize_domains_db(sample_counts_and_domains, tmp_path: Path) -> None:
    """Test visualize_domains creates visualization output image from .db database."""
    _counts_beg, _domains_beg, db_file = sample_counts_and_domains
    out_dir = tmp_path / "viz_db_out"

    visualize_domains(
        domains_path=db_file,
        counts_path=db_file,
        output_dir=out_dir,
    )

    expected_img = out_dir / "chr1_domains.png"
    assert expected_img.exists()
    assert expected_img.stat().st_size > 0
