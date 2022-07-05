import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from ..utils import _qsc_labels, info_print


__all__ = ["mini_hist", "comparison_hist"]


def mini_hist(ngenes, distribution, limit=0.02,
              ymax=1.0, mask=None, filename=None,
              x_label="Cell activation state",
              y_label="Count(%)",
              title="Distribution"):
    """
    Generates a barplot of the states with higher probabilty than the limit
    or same as the mask.
    Parameters
    ----------
    ngenes : int
        The number of genes
    distribution : np.array
        The probability distribution. Either `p_obs` or `p_out`.
    limit : float
        The threshold value to show an state.
    ymax :float
        The vertical top limit.
    mask : np.array
        An array of booleans that select which states to plot. If not provided,
        the limit would set the mask.
    filename : str
        The file name to export the image. It should have the file extension.
    x_label : str
        The label for the horizontal axis.
    y_label : str
        The label for the vertical axis.
    title : str
        The title for the barplot.
    Returns
    -------
    mask : np.array
        The selected state for a given limit. If it is already defined, the
        same array is returned
    """
    info_msg = " and exporting to {file} file.".format(file=filename) \
        if filename is not None else ""
    info_print("Plotting the {title} in a barplot{msg}"
                 .format(title=title, msg=info_msg))

    mask = distribution > limit if mask is None else mask
    length = np.sum(mask)
    labels = np.array(_qsc_labels(ngenes))
    labels = labels[mask]

    x = np.arange(length)
    y = distribution
    y = y[mask].round(3)

    plt.figure(figsize=(length / 4, 4))
    plt.bar(x, y)
    plt.xlim([-1, length])
    plt.ylim([0, ymax])
    plt.xticks(x, labels, rotation=90)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    for i, value, in enumerate(y):
        plt.text(x[i], value + 0.04 * ymax, value,
                 horizontalalignment="center",
                 verticalalignment="bottom",
                 rotation=90)

    if filename is not None:
        plt.savefig(filename, transparent=False, bbox_inches="tight")

#     plt.show()
    return mask


def comparison_hist(ngenes, p_obs, p_out, limit=0.02,
                    ymax=1.0, mask=None, filename=None,
                    x_label="Cell activation state",
                    y_label="Count(%)",
                    title="Comparison of distributions",
                    tags=["p^{obs}", "p^{out}"]):
    """
    Generates a comparison barplot between the `p_obs` and `p_out` that shows
    the states with higher probabilty than the limit or same as the mask.
    Parameters
    ----------
    ngenes : int
        The number of genes
    p_obs : np.array
        The observed probability distribution `p_obs`
    p_out : np.array
        The output probability distribution `p_out`
    limit : float
        The threshold value to show an state.
    ymax : float
        The vertical top limit.
    mask : np.array
        An array of booleans that select which states to plot. If not provided,
        the limit would set the mask.
    filenames : str
        The file name to export the image. It should have the file extension.
    x_label : str
        The label for the horizontal axis.
    y_label : str
        The label for the vertical axis.
    title : str
        The title for the barplot.
    tags : list
        Labels of each distribution in the comparison (legend).
    """
    info_msg =" and exporting to {file} file.".format(file=filename) \
        if filename is not None else ""
    info_print("Plotting the {title} in a barplot{msg}"
               .format(title=title, msg=info_msg))

    length = np.sum(mask)
    labels = np.array(_qsc_labels(ngenes))
    labels = labels[mask]

    x = 2 * np.arange(length)
    X = p_obs[mask].round(4)
    Y = p_out[mask].round(4)

    plt.figure(figsize=(length * 0.5, 4))
    plt.bar(x-0.3, X, width=0.5, label="${}$".format(tags[0]))
    plt.bar(x+0.3, Y, width=0.5, label="${}$".format(tags[1]))
    plt.xlim([-1, 2 * length])
    plt.ylim([0, ymax])
    plt.xticks(x, labels, rotation=90, fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()

    for i, (value_x, value_y) in enumerate(zip(X, Y)):
        plt.text(x[i]-0.3, value_x + 0.04*ymax, value_x,
                 horizontalalignment="center",
                 verticalalignment="bottom",
                 rotation=90, fontsize=9)

        plt.text(x[i]+0.3, value_y + 0.04*ymax, value_y,
                 horizontalalignment="center",
                 verticalalignment="bottom",
                 rotation=90, fontsize=9)

    if filename is not None:
        plt.savefig(filename, transparent=False, bbox_inches="tight")

#     plt.show()

