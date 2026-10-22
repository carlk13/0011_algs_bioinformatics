import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

# Colorblind-safe palette (Okabe-Ito inspired)
COLOR_PALETTE = {
    "primary": "#2b5c8f",  # Deep blue
    "secondary": "#c0392b",  # Crimson / Coral red
    "accent1": "#16a085",  # Teal
    "accent2": "#d35400",  # Amber / Orange
    "accent3": "#8e44ad",  # Purple
    "neutral_dark": "#2c3e50",  # Slate dark
    "neutral_light": "#ecf0f1",  # Light gray
    "gray": "#7f8c8d",  # Muted gray
}


def apply_publication_theme() -> None:
    """Configures global Matplotlib and Seaborn parameters for publication-grade figures."""
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 13,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#333333",
            "axes.linewidth": 0.8,
            "grid.color": "#e0e0e0",
            "grid.linestyle": "--",
            "grid.alpha": 0.5,
        }
    )
