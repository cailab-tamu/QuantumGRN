import os
import sys
import numpy as np
import igraph as ig
from ..qcircuit.utils import edges_to_index
from ..utils import info_print


__all__ = ["draw_network"]


def _coordinates_graph(ngenes):
    angle = 2 * np.pi / ngenes
    angles = np.arange(ngenes) * angle
    x = np.cos(angles) + 1
    y = np.sin(angles) + 1
    return x, y


def draw_network(genes, edges, theta, threshold=1, filename=None):
    """
    Draw the network representation of the QuantumGRN model.
    Parameters
    ----------
    genes : list
        The gene list from the dataset.
    edges : list
        The edges for the QuantumGRN model.
    theta : pd.Series
        The theta values in the QuantumGRN model.
    threshold : float
        The threshold for dropping edges from the network. Details in the
        manuscript
    filaname : str
        The file name to export th eimage. It should have te file extension.
    """
    info_msg = " and exporting to {file} file.".format(file=filename) \
        if filename is not None else ""
    info_print("Drawing the network representation of the qscGRN model{msg}"
                 .format(msg=info_msg))

    idx = edges_to_index(genes, edges)
    x, y = _coordinates_graph(len(genes))
    theta = theta[edges]

    theta = theta[np.abs(theta) > (threshold * np.pi / 180)]
    edges = theta.index.to_list()
    idx = edges_to_index(genes, edges)

    weigth = 30 * np.abs(theta) / np.pi
    es_color = ["#70AD47" if value >= 0 else "#FF0000" for value in theta]
    vs_color = "#DAE3F3"

    net = ig.Graph(
        n=len(genes), edges=idx,
        edge_attrs={"weigth": weigth, "color": es_color, "curved": 0.0},
        vertex_attrs={"label": genes, "color": vs_color},
        directed=False
    )

    net.vs["vertex_frame"] = "#2F528F"
    net.vs["vertex_size"] = 50
    net.vs["label_size"] = 14
    net.vs["x"] = x
    net.vs["y"] = y

    visual_style = {}
    visual_style["edge_width"] = net.es["weigth"]
    visual_style["vertex_size"] = net.vs["vertex_size"]
    visual_style["vertex_frame_color"] = net.vs["vertex_frame"]
    visual_style["vertex_shape"] = "circle"
    visual_style["bbox"] = (400, 400)
    visual_style["margin"] = 50

    ig.plot(net, filename, **visual_style)


