import matplotlib as mpl


def set_plot_style(*, column="one", font_size=8):
    """
    Configure matplotlib for publication-quality figures.

    Parameters
    ----------
    column : str
        "one" (~3.5 in) or "two" (~7.16 in) column width
    font_size : int
        Base font size in points

    Returns
    -------
    tuple
        (figure_width, figure_height) in inches
    """
    if column == "one":
        fig_w = 3.5
    elif column == "two":
        fig_w = 7.16
    else:
        raise ValueError("column must be 'one' or 'two'")

    fig_h = fig_w / 1.45

    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": font_size,
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "legend.fontsize": font_size - 1,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,

        "lines.linewidth": 1.0,
        "lines.markersize": 3.0,

        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "grid.linewidth": 0.5,
        "grid.alpha": 0.35,

        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,

        "figure.figsize": (fig_w, fig_h),
        "svg.fonttype": "none"
    })
    return (fig_w, fig_h)
