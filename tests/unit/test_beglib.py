"""Unit tests for src/beglib.py."""

from collections.abc import Iterator
from pathlib import Path

import pytest
from Bio.SeqRecord import SeqRecord

from src.beglib import read_genome


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
