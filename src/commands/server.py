"""Command handler and Dash web application for interactive sequence browser."""

import base64
import logging
import math
import sys
import tempfile
from io import StringIO
from pathlib import Path

import dash
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, callback_context, dcc, html

import src.beglib as bl

logger = logging.getLogger(__name__)


def decimate_dataframe(df: pd.DataFrame, max_points: int = 1000) -> pd.DataFrame:
    """Downsample window DataFrame if number of rows exceeds max_points.

    Args:
        df (pd.DataFrame): Input sequence window rows.
        max_points (int): Max number of window bins to return.

    Returns:
        pd.DataFrame: Downsampled DataFrame.
    """
    n = len(df)
    if n <= max_points or max_points <= 0:
        return df

    bin_size = math.ceil(n / max_points)
    group_ids = np.arange(n) // bin_size

    meta_cols = ["library_name", "seq_name", "window_start", "window_stop", "window_size"]
    count_cols = [c for c in df.columns if c not in meta_cols]

    agg_dict = {
        "library_name": "first",
        "seq_name": "first",
        "window_start": "min",
        "window_stop": "max",
        "window_size": "mean",
    }
    for col in count_cols:
        agg_dict[col] = "mean"

    decimated = df.groupby(group_ids).agg(agg_dict).reset_index(drop=True)
    return decimated


def build_figure(
    df: pd.DataFrame,
    seq_name: str,
    canonical: bool = False,
    show_n_gaps: bool = False,
    max_points: int = 1000,
) -> go.Figure:
    """Construct Plotly stacked area chart figure for selected sequence.

    Args:
        df (pd.DataFrame): Full count dataset.
        seq_name (str): Selected sequence name.
        canonical (bool): Collapse canonical k-mers if True.
        show_n_gaps (bool): Show N-gaps stack if True.
        max_points (int): Maximum bins for downsampling.

    Returns:
        go.Figure: Plotly graph object figure.
    """
    seq_df = df[df["seq_name"] == seq_name].copy()
    if seq_df.empty:
        fig = go.Figure()
        fig.update_layout(title=f"No data available for sequence: {seq_name}")
        return fig

    # Downsample raw window counts before proportion calculation for performance
    decimated_df = decimate_dataframe(seq_df, max_points=max_points)

    proc_df, feature_cols, color_map = bl.process_window_proportions(
        decimated_df,
        canonical=canonical,
        show_n_gaps=show_n_gaps,
    )

    x_vals = proc_df["window_start"].values
    fig = go.Figure()

    for col in feature_cols:
        y_vals = proc_df[col].values
        color = color_map[col]
        fig.add_trace(
            go.Scatter(
                x=x_vals,
                y=y_vals,
                mode="lines",
                name=col,
                stackgroup="one",
                line=dict(width=0.5, color=color),
                fillcolor=color,
                hovertemplate=f"<b>{col}</b>: %{{y:.3f}}<br>Pos: %{{x}}<extra></extra>",
            )
        )

    title_suffix = " (Canonical)" if canonical else ""
    if show_n_gaps:
        title_suffix += " [N-Gaps Enabled]"

    fig.update_layout(
        title=f"Genomic Signature: {seq_name}{title_suffix}",
        xaxis=dict(title="Genomic Position (bp)"),
        yaxis=dict(title="Proportion", range=[0, 1], fixedrange=True),
        template="plotly_white",
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
        ),
        margin=dict(l=60, r=160, t=60, b=60),
    )

    return fig


def create_app(
    initial_file: str | None = None,
    input_format: str = "auto",
    reloadable: bool = False,
) -> dash.Dash:
    """Create and configure the Dash browser application.

    Args:
        initial_file (str | None): Path to initial dataset file.
        input_format (str): Input format ("auto", "beg", "db").
        reloadable (bool): Allow attaching new data files in UI.

    Returns:
        dash.Dash: Configured Dash application.
    """
    app = dash.Dash(__name__, title="Birds Eye Genome - Signature Browser")

    dataset_store = {}
    if initial_file and Path(initial_file).exists():
        try:
            df = bl.load_counts_data(initial_file, file_format=input_format)
            dataset_store["file_path"] = str(Path(initial_file).resolve())
            dataset_store["sequences"] = [str(s) for s in df["seq_name"].unique()]
            dataset_store["data_json"] = df.to_json(orient="split")
        except Exception as err:
            logger.error("Failed to load initial dataset file '%s': %s", initial_file, err)

    layout_children = [
        html.Div(
            [
                html.H2("Birds Eye Genome - Sequence Signature Browser", style={"margin": "0"}),
                html.P("Interactive genomic k-mer pattern viewer", style={"color": "#666", "marginTop": "4px"}),
            ],
            style={"padding": "16px 24px", "backgroundColor": "#f8f9fa", "borderBottom": "1px solid #e9ecef"},
        ),
        dcc.Store(id="dataset-store", data=dataset_store),
    ]

    # Reload / Attachment section if reloadable
    if reloadable:
        reload_panel = html.Div(
            [
                html.Label("Attach / Load New Data File (.beg or .db):", style={"fontWeight": "bold"}),
                dcc.Upload(
                    id="upload-data",
                    children=html.Div(["Drag and Drop or ", html.A("Select BEG/DB File")]),
                    style={
                        "width": "100%",
                        "height": "50px",
                        "lineHeight": "50px",
                        "borderWidth": "1px",
                        "borderStyle": "dashed",
                        "borderRadius": "5px",
                        "textAlign": "center",
                        "marginTop": "8px",
                        "backgroundColor": "#fafafa",
                    },
                    multiple=False,
                ),
            ],
            style={"padding": "16px 24px", "borderBottom": "1px solid #eee"},
        )
        layout_children.append(reload_panel)

    controls_panel = html.Div(
        [
            html.Div(
                [
                    html.Label("Select Sequence:", style={"fontWeight": "bold", "marginBottom": "6px"}),
                    dcc.Dropdown(
                        id="seq-dropdown",
                        options=[],
                        value=None,
                        placeholder="Choose sequence...",
                        style={"width": "100%"},
                    ),
                ],
                style={"flex": "2", "minWidth": "250px", "marginRight": "16px"},
            ),
            html.Div(
                [
                    html.Label("Display Options:", style={"fontWeight": "bold", "marginBottom": "6px"}),
                    dcc.Checklist(
                        id="display-options",
                        options=[
                            {"label": " Canonical K-mers", "value": "canonical"},
                            {"label": " Show N-Gaps", "value": "n_gaps"},
                        ],
                        value=[],
                        inline=True,
                        style={"paddingTop": "6px"},
                    ),
                ],
                style={"flex": "2", "minWidth": "250px", "marginRight": "16px"},
            ),
            html.Div(
                [
                    html.Label("Downsampling (Max Windows):", style={"fontWeight": "bold", "marginBottom": "6px"}),
                    dcc.Dropdown(
                        id="decimate-dropdown",
                        options=[
                            {"label": "500 windows", "value": 500},
                            {"label": "1,000 windows (default)", "value": 1000},
                            {"label": "2,500 windows", "value": 2500},
                            {"label": "5,000 windows", "value": 5000},
                            {"label": "No downsampling (All)", "value": 0},
                        ],
                        value=1000,
                        clearable=False,
                        style={"width": "100%"},
                    ),
                ],
                style={"flex": "1", "minWidth": "200px"},
            ),
        ],
        style={
            "display": "flex",
            "flexWrap": "wrap",
            "alignItems": "center",
            "padding": "16px 24px",
            "backgroundColor": "#ffffff",
        },
    )
    layout_children.append(controls_panel)

    graph_panel = html.Div(
        [
            dcc.Loading(
                id="loading-graph",
                type="circle",
                children=dcc.Graph(id="signature-graph", style={"height": "650px"}),
            )
        ],
        style={"padding": "0 24px 24px 24px"},
    )
    layout_children.append(graph_panel)

    app.layout = html.Div(layout_children)

    # Callbacks
    if reloadable:
        @app.callback(
            Output("dataset-store", "data"),
            Input("upload-data", "contents"),
            State("upload-data", "filename"),
            prevent_initial_call=True,
        )
        def handle_upload(contents, filename):
            if not contents or not filename:
                return dash.no_update
            content_type, content_string = contents.split(",")
            decoded = base64.b64decode(content_string)

            suffix = Path(filename).suffix.lower()
            with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
                tmp.write(decoded)
                tmp_path = tmp.name

            try:
                df = bl.load_counts_data(tmp_path, file_format="auto")
                return {
                    "file_path": filename,
                    "sequences": [str(s) for s in df["seq_name"].unique()],
                    "data_json": df.to_json(orient="split"),
                }
            except Exception as err:
                logger.error("Error loading uploaded file '%s': %s", filename, err)
                return dash.no_update

    @app.callback(
        [Output("seq-dropdown", "options"), Output("seq-dropdown", "value")],
        Input("dataset-store", "data"),
    )
    def update_sequence_dropdown(store_data):
        if not store_data or "sequences" not in store_data or not store_data["sequences"]:
            return [], None
        seqs = store_data["sequences"]
        options = [{"label": seq, "value": seq} for seq in seqs]
        return options, seqs[0]

    @app.callback(
        Output("signature-graph", "figure"),
        [
            Input("seq-dropdown", "value"),
            Input("display-options", "value"),
            Input("decimate-dropdown", "value"),
        ],
        State("dataset-store", "data"),
    )
    def update_graph(selected_seq, display_opts, max_pts, store_data):
        if not selected_seq or not store_data or "data_json" not in store_data:
            fig = go.Figure()
            fig.update_layout(title="No dataset or sequence selected.")
            return fig

        df = pd.read_json(StringIO(store_data["data_json"]), orient="split")
        canonical = "canonical" in (display_opts or [])
        show_n_gaps = "n_gaps" in (display_opts or [])

        return build_figure(
            df,
            selected_seq,
            canonical=canonical,
            show_n_gaps=show_n_gaps,
            max_points=max_pts or 0,
        )

    return app


def run_server(args) -> None:
    """Execute the server command to launch the interactive browser app.

    Args:
        args: Command line arguments containing input, host, port, reloadable.
    """
    logger.info("Starting BEG interactive server...")
    input_file = getattr(args, "input", None)
    input_format = getattr(args, "input_format", "auto")
    reloadable = getattr(args, "reloadable", False)
    host = getattr(args, "host", "127.0.0.1")
    port = getattr(args, "port", 8050)

    app = create_app(
        initial_file=input_file,
        input_format=input_format,
        reloadable=reloadable,
    )

    logger.info("Serving interactive application at http://%s:%d/", host, port)
    app.run(host=host, port=port, debug=False)
