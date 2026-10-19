"""Utility script to generate hierarchical domain compositions from multi-window domain datasets."""

from __future__ import annotations

import argparse
import logging
import re
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import os

logger = logging.getLogger(__name__)


def load_domain_dataset(file_path: Path | str) -> tuple[pd.DataFrame, int]:
    """Load domain prediction data from a .beg TSV file or .db SQLite database.

    Args:
        file_path (Path | str): Path to domain file.

    Returns:
        tuple[pd.DataFrame, int]: Loaded DataFrame and inferred window size (in bp).
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Domains file not found: {file_path}")

    df = pd.DataFrame()
    if path.suffix.lower() == ".db":
        conn = sqlite3.connect(str(path))
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='domains'")
        if cursor.fetchone() is not None:
            query = """
            SELECT d.seq_name, d.window_start, d.window_stop, r.expected_states, d.domain_id, d.comment
            FROM domains d
            JOIN domain_runs r ON d.run_id = r.run_id
            ORDER BY d.seq_name, d.window_start
            """
            raw_df = pd.read_sql_query(query, conn)
            if not raw_df.empty:
                raw_df["col_name"] = "Expected_" + raw_df["expected_states"].astype(str)
                df = raw_df.pivot(
                    index=["seq_name", "window_start", "window_stop"],
                    columns="col_name",
                    values="domain_id",
                ).reset_index()
                df.columns.name = None

        # Check if windows metadata/complexity exists in SQLite
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='windows'")
        if cursor.fetchone() is not None:
            w_df = pd.read_sql_query("SELECT seq_name, window_start, window_stop, library_name, window_size, complexity FROM windows", conn)
            if not w_df.empty and not df.empty:
                df = df.merge(w_df, on=["seq_name", "window_start", "window_stop"], how="left")

        conn.close()
    else:
        df = pd.read_csv(str(path), sep="\t", comment="#")

    if df.empty:
        raise ValueError(f"Domain dataset is empty or invalid: {file_path}")

    # Determine window size
    if "window_size" in df.columns and not df["window_size"].dropna().empty:
        window_size = int(df["window_size"].dropna().iloc[0])
    elif "window_start" in df.columns and "window_stop" in df.columns:
        window_size = int((df["window_stop"] - df["window_start"] + 1).iloc[0])
    else:
        raise ValueError(f"Could not determine window size for dataset: {file_path}")

    return df, window_size


def map_parent_domain(base_df: pd.DataFrame, target_df: pd.DataFrame, col_name: str) -> list[str]:
    """Map domain labels from target_df to base_df windows based on sequence coordinate overlap.

    Args:
        base_df (pd.DataFrame): Base window DataFrame (smallest window size).
        target_df (pd.DataFrame): Target domain DataFrame.
        col_name (str): Column name in target_df (e.g. 'Expected_2').

    Returns:
        list[str]: Mapped domain label for each row in base_df.
    """
    mapped_labels = []

    # Group target_df by seq_name for fast lookup
    target_groups = {str(seq): group for seq, group in target_df.groupby("seq_name")}

    for _, row in base_df.iterrows():
        seq = str(row["seq_name"])
        b_start = row["window_start"]
        b_stop = row["window_stop"]

        if seq not in target_groups:
            mapped_labels.append("")
            continue

        seq_target = target_groups[seq]
        # Match parent window where target_start <= b_start and target_stop >= b_stop
        match = seq_target[(seq_target["window_start"] <= b_start) & (seq_target["window_stop"] >= b_stop)]
        if match.empty:
            # Fallback: calculate max coordinate overlap
            overlaps = np.maximum(
                0,
                np.minimum(b_stop, seq_target["window_stop"].values)
                - np.maximum(b_start, seq_target["window_start"].values),
            )
            best_idx = np.argmax(overlaps)
            if overlaps[best_idx] > 0:
                mapped_labels.append(str(seq_target.iloc[best_idx][col_name]))
            else:
                mapped_labels.append("")
        else:
            mapped_labels.append(str(match.iloc[0][col_name]))

    return mapped_labels


def build_hierarchical_composition(domain_files: list[str | Path]) -> pd.DataFrame:
    """Build hierarchical composition DataFrame from multi-window domain datasets.

    Args:
        domain_files (list[str | Path]): List of paths to domain output files.

    Returns:
        pd.DataFrame: Consolidated hierarchical window composition DataFrame.
    """
    if not domain_files:
        raise ValueError("At least one domain file must be provided.")

    datasets = []
    for f in domain_files:
        df, w_size = load_domain_dataset(f)
        datasets.append({"df": df, "window_size": w_size, "file": str(f)})

    # Find dataset with the smallest window size to define hierarchical base windows
    smallest_dataset = min(datasets, key=lambda x: x["window_size"])
    base_df = smallest_dataset["df"].copy().sort_values(["seq_name", "window_start"]).reset_index(drop=True)

    # Sort datasets by window size descending (Largest -> Smallest)
    sorted_datasets = sorted(datasets, key=lambda x: x["window_size"], reverse=True)

    # Extract all Expected_{N} state numbers present across all datasets
    all_states = set()
    for d in datasets:
        for c in d["df"].columns:
            if c.startswith("Expected_"):
                state_str = c.replace("Expected_", "")
                if state_str.isdigit():
                    all_states.add(int(state_str))

    sorted_states = sorted(all_states)

    # Construct result DataFrame with required BEG base metadata
    res_df = pd.DataFrame()
    res_df["library_name"] = base_df["library_name"] if "library_name" in base_df.columns else "default"
    res_df["seq_name"] = base_df["seq_name"]
    res_df["window_start"] = base_df["window_start"]
    res_df["window_stop"] = base_df["window_stop"]
    res_df["window_size"] = base_df["window_size"] if "window_size" in base_df.columns else smallest_dataset["window_size"]
    res_df["complexity"] = base_df["complexity"] if "complexity" in base_df.columns else 0.0
    res_df["entropy"] = base_df["entropy"] if "entropy" in base_df.columns else 0.0

    level_cols = []
    level_counter = 1

    # Order of levels: for each N state (ascending), iterate datasets from Largest to Smallest window size
    for n in sorted_states:
        col_name = f"Expected_{n}"
        for d in sorted_datasets:
            target_df = d["df"]
            if col_name in target_df.columns:
                level_name = f"Level_{level_counter}"
                mapped = map_parent_domain(base_df, target_df, col_name)
                res_df[level_name] = mapped
                level_cols.append(level_name)
                level_counter += 1

    # Generate shorthand hierarchical label
    shorthand_labels = []
    for _, row in res_df.iterrows():
        parts = []
        for col in level_cols:
            val = str(row[col])
            # Strip "Domain_" prefix for shorthand notation
            shorthand = val.replace("Domain_", "") if val.startswith("Domain_") else val
            parts.append(shorthand)
        shorthand_labels.append("-".join(parts))

    res_df["hierarchical_label"] = shorthand_labels

    # Reorder columns: metadata, hierarchical_label, Level_1, Level_2, ...
    final_cols = [
        "library_name",
        "seq_name",
        "window_start",
        "window_stop",
        "window_size",
        "complexity",
        "entropy",
        "hierarchical_label",
    ] + level_cols

    return res_df[final_cols]


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser for hierarchical composition tool.

    Returns:
        argparse.ArgumentParser: Configured parser.
    """
    parser = argparse.ArgumentParser(
        description="Generate hierarchical window domain composition from multi-window BEG domain outputs."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        nargs="+",
        required=True,
        help="List of path(s) to input domain prediction .beg file(s) or .db database(s).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to destination output TSV file.",
    )
    return parser


def main() -> None:
    """Main CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    try:
        res_df = build_hierarchical_composition(args.input)
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        res_df.to_csv(out_path, sep="\t", index=False)
        logger.info("Successfully generated hierarchical composition TSV with %d windows at %s", len(res_df), out_path)
    except Exception as err:
        logger.error("Failed to build hierarchical composition: %s", err)
        sys.exit(1)


if __name__ == "__main__":
    main()
