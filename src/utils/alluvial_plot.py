"""Utility script to generate alluvial plots from hierarchical composition TSV files."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def generate_alluvial_plot(
    input_file: Path | str,
    output_file: Path | str,
    title: str | None = None,
    dpi: int = 300,
) -> None:
    """Generate an alluvial plot illustrating hierarchical domain flows ordered by genomic position.

    Args:
        input_file (Path | str): Path to input hierarchical composition TSV file.
        output_file (Path | str): Destination output figure path (.png, .pdf, .svg, or .html).
        title (str, optional): Figure title string.
        dpi (int): Output figure resolution DPI (for static image formats).
    """
    in_path = Path(input_file)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    df = pd.read_csv(in_path, sep="\t", comment="#")
    if df.empty:
        raise ValueError(f"Input hierarchical composition dataset is empty: {input_file}")

    # Identify Level columns in order (Level_1, Level_2, ...)
    level_cols = [c for c in df.columns if c.startswith("Level_")]
    level_cols.sort(key=lambda x: int(x.replace("Level_", "")) if x.replace("Level_", "").isdigit() else x)

    if not level_cols:
        raise ValueError("No 'Level_{N}' columns found in hierarchical composition dataset.")

    # Sort DataFrame by seq_name and window_start to preserve genomic position order from top to bottom
    if "seq_name" in df.columns and "window_start" in df.columns:
        df = df.sort_values(["seq_name", "window_start"]).reset_index(drop=True)

    n_windows = len(df)
    n_levels = len(level_cols)

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.suffix.lower() == ".html":
        # Interactive Plotly Sankey diagram
        import plotly.graph_objects as go

        # Build node mapping per (level_idx, category_name)
        nodes = []
        node_map = {}
        node_colors = []

        # Color palette generator for domain states
        distinct_domains = set()
        for col in level_cols:
            distinct_domains.update(df[col].dropna().unique())
        sorted_domains = sorted(list(distinct_domains))
        cmap = plt.colormaps.get_cmap("tab10")
        dom_color_hex = {dom: mcolors.to_hex(cmap(i % 10)) for i, dom in enumerate(sorted_domains)}

        for l_idx, col in enumerate(level_cols):
            unique_cats = df[col].unique()
            for cat in unique_cats:
                node_id = len(nodes)
                nodes.append(f"{col}: {cat}")
                node_map[(l_idx, cat)] = node_id
                node_colors.append(dom_color_hex.get(cat, "#888888"))

        # Build links between consecutive levels
        link_source = []
        link_target = []
        link_value = []
        link_colors = []

        for l_idx in range(n_levels - 1):
            col_curr = level_cols[l_idx]
            col_next = level_cols[l_idx + 1]

            counts = df.groupby([col_curr, col_next]).size().reset_index(name="count")
            for _, row in counts.iterrows():
                c_val = row[col_curr]
                n_val = row[col_next]
                cnt = row["count"]

                src_id = node_map[(l_idx, c_val)]
                tgt_id = node_map[(l_idx + 1, n_val)]

                link_source.append(src_id)
                link_target.append(tgt_id)
                link_value.append(cnt)
                c_hex = dom_color_hex.get(c_val, "#888888")
                # Convert hex to rgba for translucent flow ribbons
                rgb = mcolors.to_rgb(c_hex)
                link_colors.append(f"rgba({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)}, 0.4)")

        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=nodes,
                color=node_colors,
            ),
            link=dict(
                source=link_source,
                target=link_target,
                value=link_value,
                color=link_colors,
            )
        )])

        fig.update_layout(
            title_text=title or "Hierarchical Domain Alluvial Plot",
            font_size=12,
            width=200 * n_levels + 300,
            height=600,
        )
        fig.write_html(str(out_path))
        logger.info("Saved interactive alluvial plot HTML to %s", out_path)

    else:
        # Static Matplotlib Alluvial Flow figure
        fig, ax = plt.subplots(figsize=(3 * n_levels + 2, max(6, n_windows * 0.15)))

        x_coords = np.arange(n_levels) * 3.0

        # Unique domain color mapping
        distinct_domains = set()
        for col in level_cols:
            distinct_domains.update(df[col].dropna().unique())
        sorted_domains = sorted(list(distinct_domains))
        cmap = plt.colormaps.get_cmap("tab10")
        dom_color_hex = {dom: cmap(i % 10) for i, dom in enumerate(sorted_domains)}

        # Height per window block (genomic position order top to bottom)
        # Y position goes from top (y = n_windows) down to 0
        y_top = np.linspace(n_windows, 0, n_windows, endpoint=False)
        window_height = (n_windows / max(1, n_windows)) * 0.85

        for i in range(n_windows):
            y_curr = y_top[i]
            row = df.iloc[i]

            for l_idx in range(n_levels):
                x = x_coords[l_idx]
                cat = row[level_cols[l_idx]]
                color = dom_color_hex.get(cat, "#CCCCCC")

                # Block for window at level l_idx
                rect = plt.Rectangle(
                    (x - 0.2, y_curr - window_height),
                    0.4,
                    window_height,
                    facecolor=color,
                    edgecolor="none",
                    alpha=0.9,
                )
                ax.add_patch(rect)

                # Connecting ribbon to next level
                if l_idx < n_levels - 1:
                    x_next = x_coords[l_idx + 1]
                    # Cubic Bezier curve flow ribbon
                    t = np.linspace(0, 1, 30)
                    x_ribbon = (1 - t) * (x + 0.2) + t * (x_next - 0.2)
                    y_ribbon_top = np.full_like(t, y_curr)
                    y_ribbon_bottom = np.full_like(t, y_curr - window_height)

                    ax.fill_between(
                        x_ribbon,
                        y_ribbon_bottom,
                        y_ribbon_top,
                        color=color,
                        alpha=0.35,
                        linewidth=0,
                    )

        ax.set_xticks(x_coords)
        ax.set_xticklabels(level_cols, fontsize=11, fontweight="bold")
        ax.set_yticks([])
        ax.set_xlim(-0.8, x_coords[-1] + 0.8)
        ax.set_ylim(-1, n_windows + 1)
        ax.set_ylabel("Genomic Positions (Top → Bottom)", fontsize=11, fontweight="bold")
        ax.set_title(title or "Hierarchical Domain Composition Alluvial Plot", fontsize=13, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        plt.tight_layout()
        plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info("Saved static alluvial plot figure to %s", out_path)


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser for alluvial plot tool.

    Returns:
        argparse.ArgumentParser: Configured parser.
    """
    parser = argparse.ArgumentParser(
        description="Generate alluvial plot from hierarchical composition TSV files."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to input hierarchical composition TSV file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to destination output figure (.png, .pdf, .svg, or .html).",
    )
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Optional figure title.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Image resolution DPI for static output (default: 300).",
    )
    return parser


def main() -> None:
    """Main CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    try:
        generate_alluvial_plot(
            input_file=args.input,
            output_file=args.output,
            title=args.title,
            dpi=args.dpi,
        )
    except Exception as err:
        logger.error("Failed to generate alluvial plot: %s", err)
        sys.exit(1)


if __name__ == "__main__":
    main()