"""Unit tests for src/utils/alluvial_plot.py."""

import pandas as pd
import pytest
from pathlib import Path

from src.utils.alluvial_plot import generate_alluvial_plot


@pytest.fixture
def sample_hierarchical_tsv(tmp_path):
    df = pd.DataFrame({
        "library_name": ["lib1", "lib1", "lib1"],
        "seq_name": ["chr1", "chr1", "chr1"],
        "window_start": [1, 1001, 2001],
        "window_stop": [1000, 2000, 3000],
        "window_size": [1000, 1000, 1000],
        "complexity": [0.5, 0.6, 0.55],
        "entropy": [2.1, 2.3, 2.2],
        "hierarchical_label": ["1-0-1-0", "1-1-1-2", "0-0-1-1"],
        "Level_1": ["Domain_1", "Domain_1", "Domain_0"],
        "Level_2": ["Domain_0", "Domain_1", "Domain_0"],
        "Level_3": ["Domain_1", "Domain_1", "Domain_1"],
        "Level_4": ["Domain_0", "Domain_2", "Domain_1"],
    })
    tsv_path = tmp_path / "hierarchical_composition.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


def test_generate_alluvial_plot_png(sample_hierarchical_tsv, tmp_path):
    out_png = tmp_path / "alluvial.png"
    generate_alluvial_plot(sample_hierarchical_tsv, out_png)
    assert out_png.exists()
    assert out_png.stat().st_size > 0


def test_generate_alluvial_plot_html(sample_hierarchical_tsv, tmp_path):
    out_html = tmp_path / "alluvial.html"
    generate_alluvial_plot(sample_hierarchical_tsv, out_html)
    assert out_html.exists()
    assert out_html.stat().st_size > 0
