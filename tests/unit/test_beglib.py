"""Unit tests for src/beglib.py."""

from collections.abc import Iterator
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.beglib import count_mers, read_genome, set_mers, window_task_generator


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
