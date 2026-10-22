"""Hybrid ReactPy web application utilizing static images and interactive Plotly tracks."""

from __future__ import annotations

import argparse
import base64
import logging
import re
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from reactpy import component, html, hooks, run

logger = logging.getLogger(__name__)

# Global cache to hold loaded dataframe and directory path
DATA_CACHE: dict = {
    "df": pd.DataFrame(),
    "img_dir": Path("."),
    "prefix": ""
}


def sanitize_filename(name: str) -> str:
    """Sanitize sequence name to match plot.py static image outputs."""
    return re.sub(r"[^\w\-.]", "_", name)


def get_image_base64(img_dir: Path, seq_name: str) -> str:
    """Load a static PNG image and encode it to Base64."""
    clean_seq = sanitize_filename(seq_name)

    # Check for direct sequence name or common output_prefix formats
    possible_paths = [
        img_dir / f"{clean_seq}.png",
        img_dir / f"{seq_name}.png"
    ]

    # If the files were generated with a prefix (e.g., 'out_sequence.png'),
    # fallback to searching for a file ending with the sequence name.
    for path in possible_paths:
        if path.exists():
            with open(path, "rb") as f:
                return base64.b64encode(f.read()).decode("utf-8")

    # Fallback search
    for file_path in img_dir.glob("*.png"):
        if clean_seq in file_path.name or seq_name in file_path.name:
            with open(file_path, "rb") as f:
                return base64.b64encode(f.read()).decode("utf-8")

    return ""


import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd


def create_interactive_tracks(seq_df: pd.DataFrame) -> go.Figure:
    """Construct a synchronized Plotly figure for metrics and stacked domain blocks."""
    if seq_df.empty:
        return go.Figure()

    level_cols = [c for c in seq_df.columns if c.startswith("Level_")]
    num_levels = len(level_cols)
    total_rows = 2 + num_levels

    # Allocate heights: Give metrics more space, divide the rest among domain tracks
    if num_levels > 0:
        row_heights = [0.1, 0.1] + [0.4 / num_levels] * num_levels
    else:
        row_heights = [0.2, 0.2]

    fig = make_subplots(
        rows=total_rows,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.03,
        row_heights=row_heights,
        subplot_titles=["Shannon Entropy", "zlib Complexity"] + [f"Hierarchy: {col}" for col in level_cols],
    )

    x_vals = seq_df["window_start"].values if "window_start" in seq_df.columns else range(len(seq_df))
    x_max = seq_df["window_stop"].max() if "window_stop" in seq_df.columns else x_vals[-1]

    # --- PANEL 1: Shannon Entropy ---
    entropy_vals = seq_df["entropy"].values if "entropy" in seq_df.columns else [0] * len(seq_df)
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=entropy_vals,
            mode="lines",
            name="Entropy",
            line={"color": "#1f77b4", "width": 1.5},
            hovertemplate="<b>Entropy</b>: %{y:.3f} nats<extra></extra>",
        ),
        row=1,
        col=1,
    )

    # --- PANEL 2: zlib Complexity ---
    complexity_vals = seq_df["complexity"].values if "complexity" in seq_df.columns else [0] * len(seq_df)
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=complexity_vals,
            mode="lines",
            name="Complexity",
            line={"color": "#d62728", "width": 1.5},
            hovertemplate="<b>Complexity</b>: %{y:.3f}<extra></extra>",
        ),
        row=2,
        col=1,
    )

    # --- PANELS 3+: Stacked Categorical Domain Blocks ---
    # Define a robust color palette for discrete domains
    domain_colors = [
        "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
        "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"
    ]

    for idx, lvl_col in enumerate(level_cols):
        row_idx = 3 + idx
        unique_domains = seq_df[lvl_col].dropna().unique()
        dom_color_map = {dom: domain_colors[i % len(domain_colors)] for i, dom in enumerate(unique_domains)}

        for _, row in seq_df.iterrows():
            w_start = row["window_start"]
            w_stop = row["window_stop"]
            dom_id = row[lvl_col]

            if pd.notna(dom_id):
                dom_id_str = str(dom_id)
                fig.add_trace(
                    go.Bar(
                        x=[w_stop - w_start],
                        y=[lvl_col],
                        base=w_start,
                        orientation="h",
                        marker_color=dom_color_map.get(dom_id, "#CCCCCC"),
                        name=dom_id_str,
                        showlegend=False,
                        hovertemplate=f"<b>Domain</b>: {dom_id_str}<br>Start: {w_start}<br>Stop: {w_stop}<extra></extra>",
                    ),
                    row=row_idx,
                    col=1,
                )

        # Clean up Y-axis for domain tracks (hide the axis text to save space)
        fig.update_yaxes(showticklabels=False, row=row_idx, col=1)

    fig.update_layout(
        height=50 + (100 * num_levels),  # Dynamically scale height based on domain tracks
        template="plotly_white",
        barmode="overlay",
        margin={"l": 60, "r": 40, "t": 40, "b": 40},
        hovermode="x unified",
        showlegend=False,
    )

    fig.update_xaxes(title_text="Genomic Position (bp)", row=total_rows, col=1, range=[0, x_max])
    fig.update_yaxes(rangemode="tozero", row=1, col=1)
    fig.update_yaxes(rangemode="tozero", row=2, col=1)

    return fig

@component
def HybridSequenceBrowser():
    """Main ReactPy UI component integrating static images and interactive plots."""
    df = DATA_CACHE["df"]
    img_dir = DATA_CACHE["img_dir"]
    prefix = DATA_CACHE["prefix"]

    sequences = sorted(df["seq_name"].astype(str).unique()) if not df.empty else []
    selected_seq, set_selected_seq = hooks.use_state(sequences[0] if sequences else "")

    if not sequences:
        return html.div({"style": {"padding": "20px"}}, "No sequence data available.")

    # Filter dataframe for selected sequence
    seq_df = df[df["seq_name"].astype(str) == selected_seq].sort_values("window_start")

    # Load Static Image
    b64_img = get_image_base64(img_dir, f"{prefix}_{selected_seq}")
    if b64_img:
        img_element = html.img(
            {
                "src": f"data:image/png;base64,{b64_img}",
                "style": {
                    "width": "100%",
                    "maxWidth": "1200px",
                    "objectFit": "contain",
                    "border": "1px solid #ddd",
                    "borderRadius": "4px",
                    "display": "block",
                },
            }
        )
    else:
        img_element = html.div(
            {
                "style": {
                    "width": "100%",
                    "maxWidth": "1200px",
                    "padding": "40px",
                    "backgroundColor": "#fee",
                    "color": "#c00",
                    "border": "1px dashed #c00",
                    "textAlign": "center",
                }
            },
            f"Static signature image for '{selected_seq}' not found in {img_dir}",
        )
    # Generate Interactive Tracks
    fig = create_interactive_tracks(seq_df)
    fig_html = fig.to_html(include_plotlyjs="cdn", full_html=False)

    return html.div(
        # ... [keep your top control bar and static image exactly as they are] ...
        html.div(
            {
                "style": {
                    "display": "flex",
                    "alignItems": "center",
                    "backgroundColor": "#ffffff",
                    "padding": "15px",
                    "borderRadius": "8px",
                    "boxShadow": "0 2px 4px rgba(0,0,0,0.1)",
                    "marginBottom": "20px",
                    "maxWidth": "1200px",
                }
            },
            html.label({"style": {"fontWeight": "bold", "marginRight": "10px"}}, "Select Sequence:"),
            html.select(
                {
                    "value": selected_seq,
                    "onChange": lambda event: set_selected_seq(event["target"]["value"]),
                    "style": {"padding": "6px 12px", "borderRadius": "4px", "minWidth": "200px"},
                },
                [html.option({"key": seq, "value": seq}, seq) for seq in sequences],
            ),
        ),
        # Static Signature Panel
        html.div({"style": {"marginBottom": "10px"}}, img_element),

        # Unified Interactive Metrics & Domain Panel
        html.iframe(
            {
                "srcDoc": fig_html,
                "style": {
                    "width": "100%",
                    "maxWidth": "1200px",
                    "height": f"{420 + (100 * len([c for c in seq_df.columns if c.startswith('Level_')]))}px",
                    "border": "none",
                    "backgroundColor": "#ffffff",
                    "borderRadius": "8px",
                    "boxShadow": "0 2px 4px rgba(0,0,0,0.1)",
                },
            }
        ),
    )


def main():
    parser = argparse.ArgumentParser(description="ReactPy hybrid server for sequence visualization.")
    parser.add_argument("-d", "--data", required=True, help="Path to TSV file containing entropy/complexity (e.g., hierarchy or counts).")
    parser.add_argument("-i", "--images", required=True, help="Directory containing pre-generated static k-mer PNG images.")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for static k-mer PNG images.")
    parser.add_argument("--host", default="127.0.0.1", help="Host address (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8000, help="Port number (default: 8000)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    img_dir = Path(args.images)
    if not img_dir.is_dir():
        logger.error("Images directory not found: %s", args.images)
        return

    logger.info("Loading sequence metrics from %s...", args.data)
    DATA_CACHE["df"] = pd.read_csv(args.data, sep="\t", comment="#")
    DATA_CACHE["img_dir"] = img_dir
    DATA_CACHE["prefix"] = args.prefix

    logger.info("Starting hybrid ReactPy server at http://%s:%d...", args.host, args.port)
    run(HybridSequenceBrowser, host=args.host, port=args.port)


if __name__ == "__main__":
    main()