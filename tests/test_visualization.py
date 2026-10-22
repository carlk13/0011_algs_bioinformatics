import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

import matplotlib.pyplot as plt

from visualization.generate_figures import (
    generate_all_figures,
    plot_figure_1_length_distribution,
    plot_figure_2_sensitivity,
    plot_figure_3_positional_error,
    plot_figure_4_error_types,
    plot_figure_5_adapter_discovery,
    plot_figure_6_demultiplexing,
)
from visualization.plot_theme import COLOR_PALETTE, apply_publication_theme


def test_plot_theme_configuration() -> None:
    """apply_publication_theme properly configures Matplotlib rcParams for publication."""
    apply_publication_theme()
    assert plt.rcParams["figure.dpi"] == 300
    assert plt.rcParams["savefig.dpi"] == 300
    assert plt.rcParams["axes.spines.top"] is False
    assert plt.rcParams["axes.spines.right"] is False
    assert "primary" in COLOR_PALETTE
    assert "secondary" in COLOR_PALETTE


def test_individual_figure_functions(tmp_path: Path) -> None:
    """Each individual plotting function produces non-empty PDF and PNG files."""
    apply_publication_theme()
    results_dir = Path("results")

    # Figure 1
    f1 = plot_figure_1_length_distribution(results_dir, tmp_path)
    assert f1.exists() and f1.stat().st_size > 0
    assert f1.with_suffix(".png").exists() and f1.with_suffix(".png").stat().st_size > 0

    # Figure 2
    f2 = plot_figure_2_sensitivity(results_dir, tmp_path)
    assert f2.exists() and f2.stat().st_size > 0
    assert f2.with_suffix(".png").exists() and f2.with_suffix(".png").stat().st_size > 0

    # Figure 3
    f3 = plot_figure_3_positional_error(results_dir, tmp_path)
    assert f3.exists() and f3.stat().st_size > 0
    assert f3.with_suffix(".png").exists() and f3.with_suffix(".png").stat().st_size > 0

    # Figure 4
    f4 = plot_figure_4_error_types(results_dir, tmp_path)
    assert f4.exists() and f4.stat().st_size > 0
    assert f4.with_suffix(".png").exists() and f4.with_suffix(".png").stat().st_size > 0

    # Figure 5
    f5 = plot_figure_5_adapter_discovery(results_dir, tmp_path)
    assert f5.exists() and f5.stat().st_size > 0
    assert f5.with_suffix(".png").exists() and f5.with_suffix(".png").stat().st_size > 0

    # Figure 6
    f6 = plot_figure_6_demultiplexing(results_dir, tmp_path)
    assert f6.exists() and f6.stat().st_size > 0
    assert f6.with_suffix(".png").exists() and f6.with_suffix(".png").stat().st_size > 0


def test_synthetic_fallback_generation(tmp_path: Path) -> None:
    """Plotting functions gracefully fall back to synthetic data when input files are missing."""
    empty_results = tmp_path / "empty_results"
    empty_results.mkdir()
    out_dir = tmp_path / "fallback_figs"

    generated = generate_all_figures(results_dir=empty_results, output_dir=out_dir)
    assert len(generated) == 6
    for path in generated.values():
        assert path.exists() and path.stat().st_size > 0
        assert (
            path.with_suffix(".png").exists()
            and path.with_suffix(".png").stat().st_size > 0
        )
