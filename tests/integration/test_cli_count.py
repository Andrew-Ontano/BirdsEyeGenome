"""Integration tests for beg count CLI command."""

import subprocess
import sys
from pathlib import Path


def test_cli_count_valid_fasta(tmp_path: Path) -> None:
    """Test running 'beg count -i dummy.fa' with a valid fasta file."""
    fasta_content = ">seq1\nATGCATGC\n>seq2\nGATTACA\n"
    fasta_file = tmp_path / "dummy.fa"
    fasta_file.write_text(fasta_content)

    result = subprocess.run(
        [sys.executable, "beg.py", "count", "-i", str(fasta_file)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0


def test_cli_count_verbose_flag(tmp_path: Path) -> None:
    """Test running 'beg count -i dummy.fa -v' outputs verbose logging to stderr."""
    fasta_content = ">seq1\nATGCATGC\n"
    fasta_file = tmp_path / "dummy.fa"
    fasta_file.write_text(fasta_content)

    result = subprocess.run(
        [sys.executable, "beg.py", "count", "-i", str(fasta_file), "-v"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "INFO" in result.stderr
    assert "Loaded accession: seq1" in result.stderr


def test_cli_count_missing_file(tmp_path: Path) -> None:
    """Test running 'beg count -i nonexistent.fa' fails cleanly with exit code 1."""
    missing_file = tmp_path / "nonexistent.fa"

    result = subprocess.run(
        [sys.executable, "beg.py", "count", "-i", str(missing_file)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "ERROR" in result.stderr
    assert "does not exist" in result.stderr


def test_cli_count_malformed_file(tmp_path: Path) -> None:
    """Test running 'beg count -i malformed.fa' fails with exit code 1 for unsupported input format."""
    malformed_file = tmp_path / "malformed.txt"
    malformed_file.write_text("This is not a fasta file format at all!")

    result = subprocess.run(
        [sys.executable, "beg.py", "count", "-i", str(malformed_file), "-f", "invalid_format_xyz"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "ERROR" in result.stderr
