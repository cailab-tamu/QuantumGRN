import numpy as np
import pandas as pd
from ..utils import info_print


__all__ = ["theta_init", "edges_init", "edges_to_index",
           "matrix_multiplication", "tensor_product"]


def theta_init(genes, method="zeros",
               activation_ratios=None):
    """
    Initialize the values of theta given a list of genes. Genes
    and activation ratios must be ordered decreasingly.
    Parameters
    ----------
    genes : list
        The list of genes from the data set.
    method : str
        The initialization method for angle parameters between
        different genes.
        zeros : Set values to 0. (Default)
        normal : Set values to a random number normally distributed.
        uniform : Set values to a random number uniformally
            distributed.
    activation_ratios : ndarray or list
        Activation values for the corresponding genes decreasingly
        ordered. The gene list must follow the same order. The
        activation ratios are use to set values in `R_y` gates on
        the encoder layer `L_enc`.
        `angle = 2 * arcsin(\sqrt(act_i))`

        If None, then `angle = \pi / 2`.
    Returns
    -------
    theta : pandas series
    """
    ngenes = len(genes)
    index = pd.MultiIndex.from_product([genes, genes],
                                       names=["control", "target"])

    if method == "zeros":
        theta = np.zeros(shape=(ngenes**2,))
    elif method == "normal":
        theta = np.random.normal(0, 1, size=(ngenes**2,))
    elif method == "uniform":
        theta = np.random.uniform(-1, 1, size=(ngenes**2,))

    sr = pd.Series(np.pi * theta, index=index)
    flag_act = "" if activation_ratios is not None else "out"

    if activation_ratios is not None:
        sqrt_activation = np.sqrt(activation_ratios)
        theta_encoder = 2 * np.arcsin(sqrt_activation)
        for i, gene in enumerate(genes):
            sr[gene, gene] = theta_encoder[i]

    else:
        for gene in genes:
            sr[gene, gene] = np.pi / 2

    info_print("Theta series is initialized using {method} as "
               "method with{flag} activation values".format(
                   method=method, flag=flag_act
               ))

    return sr


def edges_init(genes):
    """
    Creates edges for the QuantumGRN from a gene list such that
    the network is a fully-directed network. In the quantum
    circuit, it means every qubits is connected to every other
    qubit. The total number of edges is `n^2 - n`.
    Parameters
    ----------
    genes : list
        The list of genes from the data set.
    Returns
    -------
    edges : dict
    """
    combinations = []
    for g1 in genes:
        for g2 in genes:
            if g1 != g2:
                combinations.append((g1, g2))

    info_print("Edges for the QuantumGRN and quantum circuit are "
                 "created for {ngenes} genes".format(ngenes=len(genes)))

    return combinations

def edges_to_index(genes, edges):
    """
    Converts the edges dictionary (ie. edges_init) into a dictionary
    of indexes given a gene list.
    Parameters
    ----------
    genes : list
        The list of genes from the data set.
        [gene0, gene1, gene2, ...]
    edges : dict
        The list of edges with the gene name.
        [(gene0, gene1), (gene0, gene2), ...]
    Returns
    -------
    indexes = dict
        [[0, 1], [0, 2], ...]
    """
    indexes = []

    for edge in edges:
        control = genes.index(edge[0])
        target = genes.index(edge[1])
        indexes.append([control, target])

    return np.array(indexes)


def matrix_multiplication(array):
    """
    Computes the matrix multiplication for a sequences of
    transformations in a quantum circuit from left to right.
    Given a sequence T1, T2, T3, T4
    `T = T4 x T3 x T2 x T1`
    Parameters
    ----------
    array : ndaray or list
        The sequence of transformation by quantum gates. ie.
        [T1, T2, T3, T4]
    Returns
    -------
    T : ndarray
    """
    arr = None
    for matrix in array:
        arr = np.dot(matrix, arr) if arr is not None else matrix

    return arr


def tensor_product(array):
    """
    Computes the tensor multiplication for a given number of gates
    in a quantum circuit consisting of `n` qubits. The inpur array
    must follow the convention order given by Qiskit, so the result
    would match.
    `H = H4 (X) H3 (X) H2(X) H1`
    Parameters
    ----------
    array : ndarray or list
        The sequence of independent tranformation on each qubits. ie.
        [H1, H2, H3, H4]
    Returns
    -------
    H : ndarray
    """
    arr = None
    for matrix in array:
        arr = np.kron(matrix, arr) if arr is not None else matrix

    return arr

