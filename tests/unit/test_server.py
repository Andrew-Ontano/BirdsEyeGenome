"""Unit tests for src/commands/server.py."""

import pandas as pd
import plotly.graph_objects as go
from dash import Dash

from src.commands.server import build_figure, create_app, decimate_dataframe


def test_decimate_dataframe() -> None:
    """Test decimate_dataframe downsamples large window rows into bins."""
    rows = []
    for i in range(100):
        rows.append(
            {
                "library_name": "lib1",
                "seq_name": "chr1",
                "window_start": i * 10 + 1,
                "window_stop": (i + 1) * 10,
                "window_size": 10,
                "AAA": i,
                "TTT": 100 - i,
            }
        )
    df = pd.DataFrame(rows)

    # If max_points >= len(df), return unmodified df
    df_no_dec = decimate_dataframe(df, max_points=200)
    assert len(df_no_dec) == 100

    # If max_points < len(df), reduce row count
    df_dec = decimate_dataframe(df, max_points=10)
    assert len(df_dec) == 10
    assert df_dec["window_start"].iloc[0] == 1
    assert df_dec["window_stop"].iloc[-1] == 1000


def test_build_figure() -> None:
    """Test build_figure constructs Plotly go.Figure for sequence signature."""
    data = {
        "library_name": ["lib1"],
        "seq_name": ["chr1"],
        "window_start": [1],
        "window_stop": [10],
        "window_size": [10],
        "AAA": [3],
        "TTT": [3],
        "AAC": [2],
    }
    df = pd.DataFrame(data)

    fig = build_figure(df, "chr1", canonical=False, show_n_gaps=True, max_points=1000)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) > 0  # Traces present
    # Check trace names include Ns and k-mers
    trace_names = [trace.name for trace in fig.data]
    assert "Ns" in trace_names
    assert "AAA" in trace_names


def test_create_app(tmp_path) -> None:
    """Test create_app returns configured Dash application."""
    beg_file = tmp_path / "test.beg"
    lines = [
        "# Comment",
        "library_name\tseq_name\twindow_start\twindow_stop\twindow_size\tAAA\tTTT",
        "lib1\tchr1\t1\t10\t10\t5\t5",
    ]
    beg_file.write_text("\n".join(lines))

    app = create_app(initial_file=str(beg_file), reloadable=True)
    assert isinstance(app, Dash)
