"""Utility script to visualize domain predictions aligned with k-mer sequence signatures, Shannon entropy, and zlib complexity."""

from __future__ import annotations

import argparse
import logging
import re
import sqlite3
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import src.beglib as bl

logger = logging.getLogger(__name__)


def sanitize_filename(name: str) -> str:
    """Sanitize sequence name for safe use in file paths.

    Args:
        name (str): Sequence name.

    Returns:
        str: Sanitized sequence name string.
    """
    return re.sub(r"[^\w\-.]", "_", name)


def load_domain_data(file_path: Path | str) -> pd.DataFrame:
    """Load domain prediction data from a .beg TSV file or .db SQLite database.

    Args:
        file_path (Path | str): Path to domain file.

    Returns:
        pd.DataFrame: DataFrame containing `seq_name`, `window_start`, `window_stop`,
            and `Expected_{N}` state columns.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Domains file not found: {file_path}")

    if path.suffix.lower() == ".db":
        conn = sqlite3.connect(str(path))
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='domains'")
        if cursor.fetchone() is None:
            conn.close()
            return pd.DataFrame()

        query = """
        SELECT d.seq_name, d.window_start, d.window_stop, r.expected_states, d.domain_id
        FROM domains d
        JOIN domain_runs r ON d.run_id = r.run_id
        ORDER BY d.seq_name, d.window_start
        """
        raw_df = pd.read_sql_query(query, conn)
        conn.close()

        if raw_df.empty:
            return pd.DataFrame()

        raw_df["col_name"] = "Expected_" + raw_df["expected_states"].astype(str)
        pivoted = raw_df.pivot(
            index=["seq_name", "window_start", "window_stop"],
            columns="col_name",
            values="domain_id",
        ).reset_index()
        pivoted.columns.name = None
        return pivoted
    else:
        return pd.read_csv(str(path), sep="\t", comment="#")


def visualize_domains(
    domains_path: Path | str,
    counts_path: Path | str,
    output_dir: Path | str,
    canonical: bool = False,
    dpi: int = 300,
    img_format: str = "png",
) -> None:
    """Generate domain visualization figures for genomic sequences.

    Args:
        domains_path (Path | str): Path to domain predictions .beg or .db file.
        counts_path (Path | str): Path to window counts .beg or .db file.
        output_dir (Path | str): Output directory path for saving figures.
        canonical (bool): Collapse reverse complement k-mers if True.
        dpi (int): Image resolution DPI.
        img_format (str): Output image format (png, pdf, etc.).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_df = bl.load_counts_data(counts_path)
    if counts_df.empty:
        raise ValueError(f"Counts dataset is empty: {counts_path}")

    domains_df = load_domain_data(domains_path)
    if domains_df.empty:
        raise ValueError(f"Domains dataset is empty: {domains_path}")

    proc_df, feature_cols, color_map = bl.process_window_proportions(
        counts_df,
        canonical=canonical,
        show_n_gaps=False,
    )

    expected_cols = sorted(
        [c for c in domains_df.columns if c.startswith("Expected_")],
        key=lambda x: int(x.replace("Expected_", "")) if x.replace("Expected_", "").isdigit() else x,
    )

    if not expected_cols:
        raise ValueError("No 'Expected_{N}' domain prediction columns found in domains file.")

    seq_names = proc_df["seq_name"].unique()

    # Collect distinct domain IDs across all sequences to ensure consistent colors
    all_domains = set()
    for col in expected_cols:
        all_domains.update(domains_df[col].dropna().unique())
    sorted_domains = sorted([d for d in all_domains if d])

    cmap = plt.colormaps.get_cmap("tab10")
    domain_color_map = {dom: cmap(i % 10) for i, dom in enumerate(sorted_domains)}

    for seq in seq_names:
        seq_counts = proc_df[proc_df["seq_name"].astype(str) == str(seq)].sort_values("window_start")
        seq_domains = domains_df[domains_df["seq_name"].astype(str) == str(seq)].sort_values("window_start")

        if seq_counts.empty or seq_domains.empty:
            continue

        n_rows = 3 + len(expected_cols)
        fig, axes = plt.subplots(
            nrows=n_rows,
            ncols=1,
            figsize=(12, 5 + 0.8 * len(expected_cols)),
            gridspec_kw={"height_ratios": [3, 1, 1] + [0.8] * len(expected_cols)},
            sharex=True,
        )

        ax_count = axes[0]
        ax_entropy = axes[1]
        ax_complexity = axes[2]
        domain_axes = axes[3:]

        # 1. BEG stacked k-mer proportion plot
        x_vals = seq_counts["window_start"].values
        window_stops = seq_counts["window_stop"].values
        y_lower = np.zeros(len(x_vals))

        for col in feature_cols:
            y_upper = y_lower + seq_counts[col].values
            ax_count.fill_between(
                x_vals,
                y_lower,
                y_upper,
                color=color_map[col],
                label=col,
                linewidth=0.15,
            )
            y_lower = y_upper

        ax_count.set_title(f"Genomic Signature & Domain Segmentation: {seq}", fontsize=13, fontweight="bold")
        ax_count.set_ylabel("Proportion", fontsize=10)
        ax_count.set_ylim(0, 1.0)
        ax_count.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0,
            fontsize="small",
            frameon=False,
        )

        # 2. Shannon Entropy line graph
        if "entropy" in seq_domains.columns:
            entropy_vals = seq_domains["entropy"].values
        else:
            entropy_vals = bl.calculate_shannon_entropy(seq_counts[feature_cols]).values

        ax_entropy.plot(x_vals, entropy_vals, color="#1f77b4", linewidth=1.5)
        ax_entropy.set_ylabel("Entropy\n(nats)", fontsize=9, rotation=0, ha="right", va="center")
        ax_entropy.set_ylim(bottom=0)
        ax_entropy.grid(True, linestyle="--", alpha=0.5)

        # 3. zlib Complexity line graph
        if "complexity" in seq_counts.columns:
            complexity_vals = seq_counts["complexity"].values
        elif "complexity" in seq_domains.columns:
            complexity_vals = seq_domains["complexity"].values
        else:
            complexity_vals = np.zeros(len(x_vals))

        ax_complexity.plot(x_vals, complexity_vals, color="#d62728", linewidth=1.5)
        ax_complexity.set_ylabel("Complexity\n(zlib)", fontsize=9, rotation=0, ha="right", va="center")
        ax_complexity.set_ylim(bottom=0)
        ax_complexity.grid(True, linestyle="--", alpha=0.5)

        # 4. Domain state rows
        for idx, (exp_col, ax_dom) in enumerate(zip(expected_cols, domain_axes)):
            n_states = exp_col.replace("Expected_", "")
            ax_dom.set_ylabel(f"States: {n_states}", fontsize=9, rotation=0, ha="right", va="center")

            for _, row in seq_domains.iterrows():
                w_start = row["window_start"]
                w_stop = row["window_stop"]
                w_len = max(1, w_stop - w_start)
                dom_id = row.get(exp_col, "")
                if pd.notna(dom_id) and dom_id:
                    color = domain_color_map.get(dom_id, "#CCCCCC")
                    ax_dom.barh(
                        y=0,
                        width=w_len,
                        left=w_start,
                        height=0.8,
                        color=color,
                        align="center",
                        edgecolor="none",
                    )

            ax_dom.set_yticks([])
            ax_dom.set_ylim(-0.5, 0.5)

        axes[-1].set_xlabel("Genomic Position (bp)", fontsize=11)
        axes[-1].set_xlim(x_vals[0], window_stops[-1] if len(window_stops) > 0 else x_vals[0] + 1)

        # Legend for domain colors on top domain axis
        if sorted_domains and len(domain_axes) > 0:
            domain_patches = [mpatches.Patch(color=domain_color_map[dom], label=dom) for dom in sorted_domains]
            domain_axes[0].legend(
                handles=domain_patches,
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                title="Domains",
                frameon=False,
                fontsize="small",
            )

        clean_seq = sanitize_filename(str(seq))
        out_file = out_dir / f"{clean_seq}_domains.{img_format}"
        plt.tight_layout()
        plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info("Saved domain visualization figure to %s", out_file)


def build_parser() -> argparse.ArgumentParser:
    """Build command line argument parser.

    Returns:
        argparse.ArgumentParser: Argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Visualize regional domain predictions aligned with k-mer sequence signatures, entropy, and complexity."
    )
    parser.add_argument(
        "-d",
        "--domains",
        type=str,
        required=True,
        help="Path to domain predictions .beg file or .db database.",
    )
    parser.add_argument(
        "-c",
        "--counts",
        type=str,
        required=True,
        help="Path to k-mer counts .beg file or .db database.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory prefix/folder to save output visualization figures.",
    )
    parser.add_argument(
        "--canonical",
        action="store_true",
        default=False,
        help="Collapse reverse complement k-mers into canonical set.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Image resolution DPI (default: 300).",
    )
    parser.add_argument(
        "--format",
        type=str,
        default="png",
        choices=["png", "jpg", "jpeg", "pdf", "svg"],
        help="Output image format (default: png).",
    )
    return parser


def main() -> None:
    """Main execution entry point."""
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    try:
        visualize_domains(
            domains_path=args.domains,
            counts_path=args.counts,
            output_dir=args.output_dir,
            canonical=args.canonical,
            dpi=args.dpi,
            img_format=args.format,
        )
    except Exception as err:
        logger.error("Failed to visualize domains: %s", err)
        sys.exit(1)


if __name__ == "__main__":
    main()
