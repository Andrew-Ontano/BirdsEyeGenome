"""Command handler for static sequence plotting."""

import logging
import re
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

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


def run_plot(args) -> None:
    """Execute the plot command to produce static manuscript figures.

    Args:
        args: Parsed command line arguments containing input, output, and formatting settings.
    """
    logger.info("Executing 'plot' command...")
    try:
        input_path = Path(args.input)
        if not input_path.exists():
            logger.error("Input file not found: %s", args.input)
            sys.exit(1)

        input_format = getattr(args, "input_format", "auto")
        df = bl.load_counts_data(input_path, file_format=input_format)

        if df.empty:
            logger.error("Input dataset is empty: %s", args.input)
            sys.exit(1)

        output_prefix = args.output_prefix
        if not output_prefix:
            output_prefix = str(input_path.with_suffix(""))

        # Sequence filtering
        subset = getattr(args, "subset", None)
        regex = getattr(args, "regex", None)

        if subset:
            subset_list = [s.strip() for s in subset.split(",") if s.strip()]
            df = df[df["seq_name"].astype(str).isin(subset_list)]
            logger.info("Filtered dataset by subset (%d sequences match)", len(subset_list))

        if regex:
            df = df[df["seq_name"].astype(str).str.contains(regex, regex=True)]
            logger.info("Filtered dataset by regex pattern '%s'", regex)

        if df.empty:
            logger.error("No sequence records remained after filtering.")
            sys.exit(1)

        canonical = getattr(args, "canonical", False)
        show_n_gaps = getattr(args, "show_n_gaps", False)

        proc_df, feature_cols, color_map = bl.process_window_proportions(
            df,
            canonical=canonical,
            show_n_gaps=show_n_gaps,
        )

        dpi = getattr(args, "dpi", 300)
        output_format = getattr(args, "output_format", "png")
        key_setting = getattr(args, "key", "yes").lower()

        seq_names = proc_df["seq_name"].unique()
        logger.info("Generating static plots for %d sequence(s)...", len(seq_names))

        for seq in seq_names:
            seq_df = proc_df[proc_df["seq_name"] == seq]
            x_vals = seq_df["window_start"].values

            fig, ax = plt.subplots(figsize=(12, 5))

            y_lower = np.zeros(len(x_vals))
            for col in feature_cols:
                y_upper = y_lower + seq_df[col].values
                ax.fill_between(
                    x_vals,
                    y_lower,
                    y_upper,
                    color=color_map[col],
                    label=col,
                    linewidth=.15,
                    edgecolor=color_map[col]
                )
                y_lower = y_upper

            ax.set_title(f"Genomic Signature: {seq}", fontsize=14, fontweight="bold")
            ax.set_xlabel("Genomic Position (bp)", fontsize=11)
            ax.set_ylabel("Proportion", fontsize=11)
            ax.set_ylim(0, 1.0)
            ax.set_xlim(x_vals[0], x_vals[-1] if len(x_vals) > 1 else x_vals[0] + 1)

            if key_setting == "yes":
                ax.legend(
                    bbox_to_anchor=(1.02, 1),
                    loc="upper left",
                    borderaxespad=0,
                    fontsize="small",
                    frameon=False,
                )

            clean_seq = sanitize_filename(str(seq))
            out_filename = f"{output_prefix}_{clean_seq}.{output_format}"
            plt.savefig(out_filename, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            logger.info("Saved static plot to %s", out_filename)

        if key_setting == "alone":
            fig_key, ax_key = plt.subplots(figsize=(8, 6))
            ax_key.axis("off")
            patches = [mpatches.Patch(color=color_map[col], label=col) for col in feature_cols]
            ax_key.legend(
                handles=patches,
                loc="center",
                ncol=min(8, max(1, len(feature_cols) // 8 + 1)),
                frameon=False,
                fontsize="medium",
            )
            key_filename = f"{output_prefix}_key.{output_format}"
            plt.savefig(key_filename, dpi=dpi, bbox_inches="tight")
            plt.close(fig_key)
            logger.info("Saved standalone key image to %s", key_filename)

        logger.info("Plot command completed successfully.")

    except Exception as err:
        logger.error("Failed to generate plots: %s", err)
        sys.exit(1)
