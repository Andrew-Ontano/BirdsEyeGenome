"""Unit tests for src/beglib.py."""

from collections.abc import Iterator
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.beglib import (
    calculate_complexity,
    count_mers,
    generate_kmer_colors,
    get_canonical_kmer,
    get_reverse_complement,
    init_db,
    load_counts_data,
    process_window_proportions,
    read_genome,
    set_mers,
    window_task,
    window_task_generator,
    write_counts_beg,
    write_counts_db,
)


def test_read_genome_valid_fasta(tmp_path: Path) -> None:
    """Test reading a valid FASTA file returns SeqRecord iterator."""
    fasta_content = ">seq1 First chromosome\nATGCATGC\n>seq2 Second chromosome\nCGTACS\n"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    records = read_genome(fasta_file, input_format="fasta")
    assert isinstance(records, Iterator)

    record_list = list(records)
    assert len(record_list) == 2
    assert isinstance(record_list[0], SeqRecord)
    assert record_list[0].id == "seq1"
    assert str(record_list[0].seq) == "ATGCATGC"
    assert record_list[1].id == "seq2"


def test_read_genome_file_not_found(tmp_path: Path) -> None:
    """Test that read_genome raises FileNotFoundError for missing files."""
    missing_file = tmp_path / "nonexistent.fasta"
    with pytest.raises(FileNotFoundError, match="does not exist or is not a file"):
        read_genome(missing_file)


def test_read_genome_directory(tmp_path: Path) -> None:
    """Test that read_genome raises FileNotFoundError if path is a directory."""
    with pytest.raises(FileNotFoundError, match="does not exist or is not a file"):
        read_genome(tmp_path)


def test_set_mers_default() -> None:
    """Test set_mers generates 64 canonical trinucleotides."""
    kmers = set_mers(3)
    assert len(kmers) == 64
    assert kmers[0] == "AAA"
    assert kmers[-1] == "TTT"


def test_set_mers_custom_kmer_size() -> None:
    """Test set_mers with custom kmer size."""
    kmers_2 = set_mers(2)
    assert len(kmers_2) == 16
    assert kmers_2 == sorted(kmers_2)


def test_window_task_generator() -> None:
    """Test window task generator correctly creates sliding window tasks."""
    records = [
        SeqRecord(Seq("ATGCATGCAT"), id="seq1"),
    ]
    kmers = set_mers(3)
    # seq_len=10, window_size=4, window_interval=3 -> starts at 0, 3, 6, 9
    tasks = list(window_task_generator(records, window_size=4, window_interval=3, kmer_list=kmers))
    assert len(tasks) == 4

    assert tasks[0]["seq_name"] == "seq1"
    assert tasks[0]["start"] == 0
    assert tasks[0]["end"] == 4
    assert tasks[0]["window_size"] == 4
    assert tasks[0]["sequence"] == "ATGC"

    assert tasks[2]["start"] == 6
    assert tasks[2]["end"] == 10
    assert tasks[2]["window_size"] == 4
    assert tasks[2]["sequence"] == "GCAT"

    assert tasks[3]["start"] == 9
    assert tasks[3]["end"] == 10
    assert tasks[3]["window_size"] == 1
    assert tasks[3]["sequence"] == "T"


def test_calculate_complexity() -> None:
    """Test calculate_complexity returns valid zlib compression ratio."""
    assert calculate_complexity("") == 0.0
    comp_poly = calculate_complexity("A" * 100)
    comp_rand = calculate_complexity("ATGC" * 25)
    assert comp_poly > 0.0
    assert comp_rand > 0.0


def test_count_mers() -> None:
    """Test count_mers correctly tallies overlapping k-mers in a window."""
    kmers = set_mers(2)
    task = {
        "seq_name": "seq1",
        "start": 0,
        "end": 4,
        "window_size": 4,
        "sequence": "AAAA",
        "kmers": kmers,
    }
    result = count_mers(task)
    assert result["seq_name"] == "seq1"
    assert result["counts"]["AA"] == 3
    assert result["counts"]["AC"] == 0


def test_window_task() -> None:
    """Test window_task returns counts and complexity."""
    kmers = set_mers(2)
    task = {
        "seq_name": "seq1",
        "start": 0,
        "end": 4,
        "window_size": 4,
        "sequence": "AAAA",
        "kmers": kmers,
    }
    result = window_task(task)
    assert result["seq_name"] == "seq1"
    assert "complexity" in result
    assert result["complexity"] > 0.0
    assert result["counts"]["AA"] == 3


def test_init_db_and_write_counts_db(tmp_path: Path) -> None:
    """Test init_db creates table and write_counts_db populates sqlite table."""
    db_file = tmp_path / "output.db"
    kmers = ["AA", "AC", "AG", "AT"]
    conn = init_db(db_file, kmers)

    results = [
        {
            "seq_name": "seq1",
            "start": 0,
            "end": 10,
            "window_size": 10,
            "counts": {"AA": 2, "AC": 1, "AG": 0, "AT": 0},
        }
    ]
    write_counts_db(conn, "test.fa", results, kmers)

    cursor = conn.cursor()
    cursor.execute("SELECT library_name, seq_name, window_start, window_stop, window_size FROM windows")
    rows = cursor.fetchall()
    assert len(rows) == 1
    assert rows[0] == ("test.fa", "seq1", 1, 10, 10)

    cursor.execute("SELECT key, value FROM metadata")
    metadata = dict(cursor.fetchall())
    assert metadata["kmer_size"] == "2"
    assert metadata["library_name"] == "test.fa"

    conn.close()

    df = load_counts_data(db_file)
    assert len(df) == 1
    assert df["AA"].iloc[0] == 2
    assert df["AC"].iloc[0] == 1


def test_sparse_storage_and_kmer_6(tmp_path: Path) -> None:
    """Test normalized db structure with k-mer size 6 sparse counting."""
    kmers_6 = set_mers(6)
    db_file = tmp_path / "kmer6.db"
    conn = init_db(db_file, kmers_6, kmer_size=6)

    results = [
        {
            "seq_name": "chr1",
            "start": 0,
            "end": 100,
            "window_size": 100,
            "counts": {"AAAAAA": 5, "CCCCCC": 2},
        }
    ]
    write_counts_db(conn, "genome6.fa", results, kmers_6, kmer_size=6)
    conn.close()

    df = load_counts_data(db_file)
    assert len(df) == 1
    assert df["library_name"].iloc[0] == "genome6.fa"
    assert df["seq_name"].iloc[0] == "chr1"
    assert df["AAAAAA"].iloc[0] == 5
    assert df["CCCCCC"].iloc[0] == 2
    assert df["GGGGGG"].iloc[0] == 0
    assert len(df.columns) == 6 + 4096


def test_write_counts_beg(tmp_path: Path) -> None:
    """Test write_counts_beg outputs TSV with header comment metadata."""
    beg_file = tmp_path / "output.beg"
    kmers = ["AA", "AC"]
    results = [
        {
            "seq_name": "seq1",
            "start": 0,
            "end": 10,
            "window_size": 10,
            "counts": {"AA": 2, "AC": 1},
        }
    ]
    write_counts_beg(beg_file, tmp_path / "dummy.fa", "dummy.fa", kmers, results)

    lines = beg_file.read_text().splitlines()
    assert lines[0].startswith("# Input file:")
    assert lines[1].startswith("# Generated by Birds Eye Genome")
    assert lines[2] == "library_name\tseq_name\twindow_start\twindow_stop\twindow_size\tcomplexity\tAA\tAC"
    assert lines[3] == "dummy.fa\tseq1\t1\t10\t10\t0.0\t2\t1"


def test_generate_kmer_colors() -> None:
    """Test generate_kmer_colors produces valid hex colors for each k-mer."""
    kmers = set_mers(2)
    colors = generate_kmer_colors(kmers)
    assert len(colors) == len(kmers)
    for hex_code in colors.values():
        assert hex_code.startswith("#")
        assert len(hex_code) == 7


def test_canonical_and_reverse_complement() -> None:
    """Test reverse complement and canonical k-mer helper functions."""
    assert get_reverse_complement("AAA") == "TTT"
    assert get_reverse_complement("ATC") == "GAT"
    assert get_canonical_kmer("TTT") == "AAA"
    assert get_canonical_kmer("ATC") == "ATC"
    assert get_canonical_kmer("GAT") == "ATC"


def test_load_counts_data_beg_and_db(tmp_path: Path) -> None:
    """Test loading window counts from .beg and .db files."""
    kmers = ["AAA", "AAC"]
    results = [
        {
            "seq_name": "chr1",
            "start": 0,
            "end": 10,
            "window_size": 10,
            "counts": {"AAA": 3, "AAC": 1},
        }
    ]

    beg_file = tmp_path / "test.beg"
    write_counts_beg(beg_file, "dummy.fa", "dummy.fa", kmers, results)

    df_beg = load_counts_data(beg_file)
    assert len(df_beg) == 1
    assert df_beg["seq_name"].iloc[0] == "chr1"
    assert df_beg["AAA"].iloc[0] == 3

    db_file = tmp_path / "test.db"
    conn = init_db(db_file, kmers)
    write_counts_db(conn, "dummy.fa", results, kmers)
    conn.close()

    df_db = load_counts_data(db_file)
    assert len(df_db) == 1
    assert df_db["seq_name"].iloc[0] == "chr1"
    assert df_db["AAA"].iloc[0] == 3


def test_process_window_proportions() -> None:
    """Test process_window_proportions with canonical mapping and N-gaps."""
    data = {
        "library_name": ["lib1"],
        "seq_name": ["chr1"],
        "window_start": [1],
        "window_stop": [10],
        "window_size": [10],
        "AAA": [2],
        "TTT": [2],
        "AAC": [1],
        "GTT": [1],
    }
    df = pd.DataFrame(data)

    # Standard proportions without N-gaps or canonical collapse
    df_proc, features, _color_map = process_window_proportions(df, canonical=False, show_n_gaps=False)
    assert set(features) == {"AAA", "TTT", "AAC", "GTT"}
    assert pytest.approx(df_proc["AAA"].iloc[0]) == 2 / 6
    assert pytest.approx(df_proc["TTT"].iloc[0]) == 2 / 6

    # Canonical collapse
    df_canon, canon_features, _ = process_window_proportions(df, canonical=True, show_n_gaps=False)
    assert "TTT" not in canon_features
    assert "AAA" in canon_features
    # AAA + TTT = 4 / 6
    assert pytest.approx(df_canon["AAA"].iloc[0]) == 4 / 6

    # N-gaps calculation: Ns = max(0, window_size - sum(kmers) - kmer_len + 1)
    # window_size=10, sum(kmers)=6, kmer_len=3 -> Ns = 10 - 6 - 3 + 1 = 2
    # total = 6 + 2 = 8
    df_ngaps, ngap_features, ngap_colors = process_window_proportions(df, canonical=False, show_n_gaps=True)
    assert ngap_features[0] == "Ns"
    assert ngap_colors["Ns"] == "#000000"
    assert pytest.approx(df_ngaps["Ns"].iloc[0]) == 2 / 8
    assert pytest.approx(df_ngaps["AAA"].iloc[0]) == 2 / 8
