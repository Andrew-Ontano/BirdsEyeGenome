"""Birds Eye Genome core library module."""

import itertools
import sqlite3
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_genome(
    input_file: str | Path,
    input_format: str = "fasta",
) -> Iterator[SeqRecord]:
    """Check that the input file exists and load sequence records via Biopython SeqIO.

    Args:
        input_file (str | Path): Path to the input genomic sequence file.
        input_format (str, optional): Input sequence file format recognized by SeqIO.
            Defaults to "fasta".

    Returns:
        Iterator[SeqRecord]: An iterator yielding SeqRecord objects from the input file.

    Raises:
        FileNotFoundError: If the specified file does not exist or is not a regular file.
    """
    file_path = Path(input_file)
    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(f"Input file does not exist or is not a file: {input_file}")

    return SeqIO.parse(str(file_path), input_format)


def set_mers(kmer_size: int = 3) -> list[str]:
    """Build all possible k-mers constructed from A, C, G, and T.

    Args:
        kmer_size (int, optional): Length of the k-mers to generate. Defaults to 3.

    Returns:
        list[str]: Sorted list of all k-mer strings constructed from A, C, G, T.
    """
    nucleotides = ["A", "C", "G", "T"]
    kmers = ["".join(p) for p in itertools.product(nucleotides, repeat=kmer_size)]
    return sorted(kmers)


def window_task_generator(
    records: Iterable[SeqRecord],
    window_size: int,
    window_interval: int,
    kmer_list: list[str],
) -> Iterator[dict[str, Any]]:
    """Define sliding windows for each genomic sequence and yield window tasks.

    Args:
        records (Iterable[SeqRecord]): Iterable yielding SeqRecord objects.
        window_size (int): The window size in base pairs.
        window_interval (int): The interval step size between sliding window start positions.
        kmer_list (list[str]): The list of k-mers to be counted.

    Yields:
        dict[str, Any]: A dictionary containing window task parameters:
            - "seq_name": Sequence accession/ID.
            - "start": Start position (0-based inclusive).
            - "end": End position (0-based exclusive).
            - "window_size": Actual window length in base pairs.
            - "sequence": Subsequence string for the window.
            - "kmers": List of k-mer strings to count.
    """
    for record in records:
        seq_str = str(record.seq).upper()
        seq_len = len(seq_str)
        if seq_len == 0:
            continue

        for start in range(0, seq_len, window_interval):
            end = min(start + window_size, seq_len)
            subseq = seq_str[start:end]
            actual_window_size = len(subseq)
            yield {
                "seq_name": str(record.id),
                "start": start,
                "end": end,
                "window_size": actual_window_size,
                "sequence": subseq,
                "kmers": kmer_list,
            }


def count_mers(task: dict[str, Any]) -> dict[str, Any]:
    """Count occurrences of each k-mer within a window sequence task.

    Args:
        task (dict[str, Any]): Window task dictionary containing "seq_name", "start",
            "end", "window_size", "sequence", and "kmers".

    Returns:
        dict[str, Any]: Dictionary containing task metadata and a "counts" mapping of
            each k-mer to its count in the window.
    """
    subseq = task["sequence"].upper()
    kmers = task["kmers"]
    kmer_len = len(kmers[0]) if kmers else 0

    counts = {kmer: 0 for kmer in kmers}
    seq_len = len(subseq)

    if kmer_len > 0 and seq_len >= kmer_len:
        for i in range(seq_len - kmer_len + 1):
            kmer = subseq[i : i + kmer_len]
            if kmer in counts:
                counts[kmer] += 1

    return {
        "seq_name": task["seq_name"],
        "start": task["start"],
        "end": task["end"],
        "window_size": task["window_size"],
        "counts": counts,
    }


def init_db(db_path: str | Path, kmers: list[str]) -> sqlite3.Connection:
    """Initialize SQLite database with a tabulation table for k-mer counts.

    Args:
        db_path (str | Path): Path to SQLite database file.
        kmers (list[str]): List of k-mer names to include as count columns.

    Returns:
        sqlite3.Connection: Active connection to the created database.
    """
    db_file = Path(db_path)
    if db_file.exists():
        db_file.unlink()

    conn = sqlite3.connect(str(db_file))
    cursor = conn.cursor()

    kmer_cols = ", ".join([f'"{kmer}" INTEGER' for kmer in kmers])
    create_table_sql = f"""
    CREATE TABLE IF NOT EXISTS kmer_counts (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        library_name TEXT,
        seq_name TEXT,
        window_start INTEGER,
        window_stop INTEGER,
        window_size INTEGER
        {", " + kmer_cols if kmer_cols else ""}
    )
    """
    cursor.execute(create_table_sql)
    conn.commit()
    return conn


def write_counts_db(
    conn: sqlite3.Connection,
    library_name: str,
    results: list[dict[str, Any]],
    kmers: list[str],
) -> None:
    """Insert window count results into the SQLite database table.

    Args:
        conn (sqlite3.Connection): Active SQLite database connection.
        library_name (str): Shortened library name (e.g. filename).
        results (list[dict[str, Any]]): List of window count result dictionaries.
        kmers (list[str]): List of k-mers corresponding to count columns.
    """
    cursor = conn.cursor()
    columns = ["library_name", "seq_name", "window_start", "window_stop", "window_size"] + [f'"{k}"' for k in kmers]
    placeholders = ", ".join(["?"] * len(columns))
    col_str = ", ".join(columns)
    sql = f"INSERT INTO kmer_counts ({col_str}) VALUES ({placeholders})"

    rows = []
    for res in results:
        counts = res.get("counts", {})
        row = [
            library_name,
            res["seq_name"],
            res["start"],
            res["end"],
            res["window_size"],
        ] + [counts.get(kmer, 0) for kmer in kmers]
        rows.append(row)

    cursor.executemany(sql, rows)
    conn.commit()


def write_counts_beg(
    output_file: str | Path,
    input_file: str | Path,
    library_name: str,
    kmers: list[str],
    results: list[dict[str, Any]],
) -> None:
    """Write window count results to a tab-separated BEG file with comment metadata.

    Args:
        output_file (str | Path): Path to destination BEG file.
        input_file (str | Path): Total path to the input sequence file.
        library_name (str): Shortened library name.
        kmers (list[str]): List of k-mers.
        results (list[dict[str, Any]]): List of window count result dictionaries.
    """
    abs_input_path = Path(input_file).resolve()
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(f"# Input file: {abs_input_path}\n")
        f.write("# Generated by Birds Eye Genome (beg count)\n")

        header = ["library_name", "seq_name", "window_start", "window_stop", "window_size"] + kmers
        f.write("\t".join(header) + "\n")

        for res in results:
            counts = res.get("counts", {})
            row = [
                library_name,
                res["seq_name"],
                str(res["start"]),
                str(res["end"]),
                str(res["window_size"]),
            ] + [str(counts.get(kmer, 0)) for kmer in kmers]
            f.write("\t".join(row) + "\n")
