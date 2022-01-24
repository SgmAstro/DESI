'''
Generate publication ready plots
'''
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def set_size(width=240, fraction=1, equal=False):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Parameters
    ----------
    width: float
            Width in pts
    fraction: float
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    if equal == False:
        golden_ratio = (5 ** 0.5 - 1) / 2
    else:
        # height = width
        golden_ratio = 1.

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio

    return (fig_width_in, fig_height_in)


plt.style.use(os.environ['CODE_ROOT'] + '/style.mplstyle')

nice_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "font.serif" : "Times New Roman",
}

matplotlib.rcParams['figure.figsize'] = set_size()
matplotlib.rcParams.update(nice_fonts)

