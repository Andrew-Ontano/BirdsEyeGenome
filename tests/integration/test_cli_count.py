"""Integration tests for beg count CLI command."""

import sqlite3
import subprocess
import sys
from pathlib import Path


def test_cli_count_beg_output(tmp_path: Path) -> None:
    """Test running 'beg count' creates a valid BEG TSV output file."""
    fasta_content = ">seq1\nATGCATGC\n"
    fasta_file = tmp_path / "dummy.fa"
    fasta_file.write_text(fasta_content)
    output_beg = tmp_path / "out.beg"

    result = subprocess.run(
        [sys.executable, "beg.py", "count", "-i", str(fasta_file), "-o", str(output_beg)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert output_beg.exists()
    content = output_beg.read_text()
    assert "# Input file:" in content
    assert "library_name\tseq_name\twindow_start\twindow_stop\twindow_size" in content


def test_cli_count_db_output(tmp_path: Path) -> None:
    """Test running 'beg count -p db' creates a valid SQLite database output."""
    fasta_content = ">seq1\nATGCATGC\n"
    fasta_file = tmp_path / "dummy.fa"
    fasta_file.write_text(fasta_content)
    output_db = tmp_path / "out.db"

    result = subprocess.run(
        [sys.executable, "beg.py", "count", "-i", str(fasta_file), "-o", str(output_db), "-p", "db"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert output_db.exists()

    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()
    cursor.execute("SELECT library_name, seq_name FROM kmer_counts")
    rows = cursor.fetchall()
    assert len(rows) == 1
    assert rows[0] == ("dummy.fa", "seq1")
    conn.close()


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
    assert "Successfully completed k-mer counting" in result.stderr


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
