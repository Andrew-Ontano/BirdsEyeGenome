"""Birds Eye Genome core library module."""

import itertools
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
