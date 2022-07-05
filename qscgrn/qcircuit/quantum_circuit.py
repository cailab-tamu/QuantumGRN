import numpy as np
import pandas as pd
from .utils import *
from .gates import *
from ..run import *
from ..utils import info_print


__all__ = ['quantum_circuit']


class quantum_circuit(qscgrn_model):
    """
    Attributes
    ----------
    ngenes : int
        Number of genes in the Quanutm GRN.
    genes : list
        Gene list for QuantumGRN modelling.
    theta : pd.Series
        Theta values given edges in the QuantumGRN.
    edges : list
        Edges for the QuantumGRN.
    indexes : list
        Numerical values of the edges for usage in the quantum circuit.
    encoder : np.array
        Matrix representation of the encoder layer `L_enc`.
    regulation : np.array
        Matrix representation of the regulation layers `L_k`. The
        encoder layers are grouped in a single array.
    circuit : np.array
        A single matrix representation for the quantum circuit
        transformation.
    input : np.array
        The quantum state as input for the quantum circuit model.
    derivatives : pd.DataFrame
        Derivatives of the quantum state in the output register
        with respect of each parameter.
    drop_zero : bool
        If True, a normalization step for `p^out` that set the `|0>_n`
        state to 0, and rescale the rest of the distribution.

    Methods
    -------
    computer_encoder()
        Computes the transformation matrix for the `L_enc` into the
        encoder attribute.
    compute_regulation()
        Computes the transformation matrix for the `L_k` into the
        regulation attribute.
    generate_circuit()
        Compute the transformation matrix for the `L_enc` and `L_k`.
    tranform_matrix()
        Compute the transformation matrix for the entire quantum
        circuit.
    output_state()
        Compute the quantum state in the output register given an
        input state.
    output_probabilities(drop_zeros)
        Compute the probability distribution in the output register.
        If drop_zero is True, a normalization step is done.
    create_derivatives()
        Creates a pd.DataFrame to store the derivatives of the
        output state with respect of the parameters.
    der_encoder()
        Computes the derivatives with respect of the parameters
        in the `L_enc` layer.
    der_regulation()
        Computes the derivatives with respect of the parameters
        in the `L_k` layers.
    compute_derivates()
        Computes the derivatives by calling the der_encoder
        and the der_regulation methods.
    """

    def __init__(self, genes, theta, edges, drop_zero=True):
        """
        Parameters
        ----------
        genes : list
            Gene list for QuantumGRN modelling.
        theta : pd.Series
            Theta values given edges in the QuantumGRN.
        edges : list
            Edges for the QuantumGRN.
        drop_zero : bool
            If True, a normalization step for `p^out` that set the
            `|0>_n` state to 0, and rescale the rest of the
            distribution.
        """
        super().__init__(genes, theta, edges, drop_zero)
        # numerical indexes are needed to construct the circuit
        self.indexes = edges_to_index(genes, edges)
        # array storage for the quantum circuit (Lenc, Lk and
        # and transformation matrix)
        self.encoder = None
        self.regulation = None
        self.circuit = False
        # parameters for quantum circuit such as input state,
        # derivatives and drop_zero
        self.input = np.zeros((2**self.ngenes, 1))
        self.input[0, 0] = 1.
        self.derivatives = None


    def __str__(self):
        return ("Quantum circuit for {ngenes} genes for GRN"
                "modelling".format(ngenes=len(self.genes)))


    def _circuit_is_empty(self):
        """
        Validates whether the quantum circuit is initialized or not.
        Raises
        ------
        AttributeError
            If circuit attribute is a None object.
        """
        if not self.circuit:
            info_print("The QuantumGRN model is not initialized")
            raise AttributeError("The quantum circuit for GRN model "
                                 "is not constructed")


    def _der_is_not_empty(self):
        """
        Validates if the derivatives for the quantum circuit are
        not initialized.
        Raises
        ------
        AttributeError
            If derivatives is not a None object
        """
        if self.derivatives is not None:
            info_print("Derivatives for the QuantumGRN are already "
                          "initialized", level="E")
            raise AttributeError("The quantum circuit for GRN model "
                                  "has derivatives initialized")


    def _der_is_empty(self):
        """
        Validates if the derivatives for the quantum circuit are
        initialized
        Raises
        ------
        AttributeError
            If derivatives is a None object
        """
        if self.derivatives is None:
            info_print("Derivatives for the QuantumGRN are not "
                          "initialized", level="E")
            raise AttributeError("The quantum circuit for GRN model "
                                 "does not have derivatives "
                                 "initialized")


    def compute_encoder(self):
        """
        Computes the transformation matrices of each gate in `L_enc`
        layer and saves the result into self.encoder
        """
        RR = np.zeros((len(self.genes), 2, 2))

        for idx, gene in enumerate(self.genes):
            RR[idx] = ry_gate(self.theta[(gene, gene)])

        self.encoder = RR


    def compute_regulation(self):
        """
        Computes the tranformation matrices of each gate in `L_k`
        layer and saves the result into self.regulation
        """
        arr = np.zeros((len(self.edges), \
                        2**self.ngenes, 2**self.ngenes))

        for i, edge in enumerate(self.edges):
            idx = self.indexes[i]
            arr[i] = cry_gate(self.theta[edge], self.ngenes,
                              idx[0], idx[1])

        self.regulation = arr


    def generate_circuit(self):
        """
        Computes the `L_enc` and `L_k` accordingly to parameters
        such as `theta` and edges.
        """
        self.circuit = False
        self.compute_encoder()
        self.compute_regulation()
        self.circuit = True


    def transform_matrix(self):
        """
        Computes the tranformation matrix for the quantum circuit
        once `L_enc` and `L_k` are computed.
        Returns
        -------
        T : ndarray
            The transformation matrix for the whole quantum circuit.
        """
        self._circuit_is_empty()
        transform_regulation = matrix_multiplication(self.regulation)
        transform_encoder = tensor_product(self.encoder)
        return np.dot(transform_regulation, transform_encoder)


    def output_state(self):
        """
        Computes the quantum state in the output register of the
        quantum circuit given an input state.
        Returns
        -------
        `psi_out` : ndarray
            The quantum state in the output register
        """
        self._circuit_is_empty()
        return np.dot(self.transform_matrix(), self.input)


    def output_probabilities(self, drop_zero):
        """
        Computes the prbability distribution in the output register
        of the quantum circuit given an input state.
        Parameters
        ----------
        drop_zero : bool
            If True, the probability of `|0>_n` is set to 0, and the
            rest of the distribution is rescaled.
        Returns
        -------
        `p^out` : ndarray
            The output probability distribution `p^out` whether it
            was normalized or not.
        """
        self._circuit_is_empty()
        probabilities = np.power(self.output_state(), 2)

        if drop_zero:
            probabilities[0, 0] = 0.
            return probabilities / np.sum(probabilities)

        else:
            return probabilities


    def create_derivatives(self):
        """
        Creates the derivatives of the quantum state in the output
        register with respect of each parameter.
        """
        self._der_is_not_empty()
        index = pd.MultiIndex.from_product([self.genes, self.genes],
                                           names=["control", "target"])
        self.derivatives = pd.DataFrame(np.zeros((self.ngenes**2,
                                                  2**self.ngenes)),
                                        index=index)


    def der_encoder(self):
        """
        Computes the derivatives with respect of the parameters
        in the `L_enc` layer and save it in the derivative
        attribute.
        """
        self._circuit_is_empty()
        self._der_is_empty()

        for idx, gene in enumerate(self.genes):
            der_matrix = np.array(self.encoder)
            gate = der_ry_gate(self.theta[(gene, gene)])
            der_matrix[idx] = gate

            transform_regulation = matrix_multiplication(self.regulation)
            transform_encoder = tensor_product(der_matrix)
            transform = np.dot(transform_regulation, transform_encoder)

            derivative = np.dot(transform, self.input)
            self.derivatives.loc[(gene, gene)] = \
                    derivative.reshape(2**self.ngenes,)


    def der_regulation(self):
        """
        Computes the derivatives with respect of the parameters
        in the `L_k` layers and save it in the derivative
        attribute.
        """
        self._circuit_is_empty()
        self._der_is_empty()

        for idx, edge in enumerate(self.edges):
            index = self.indexes[idx]
            der_matrix = np.array(self.regulation)
            gate = der_cry_gate(self.theta[edge], self.ngenes,
                                index[0], index[1])
            der_matrix[idx] = gate

            transform_regulation = matrix_multiplication(der_matrix)
            transform_encoder = tensor_product(self.encoder)
            transform = np.dot(transform_regulation, transform_encoder)

            derivative = np.dot(transform, self.input)
            self.derivatives.loc[edge] = \
                derivative.reshape(2**self.ngenes,)


    def compute_derivatives(self):
        """
        Computes the derivatives by calling the der_encoder and
        der_regulation methods.
        """
        self.derivatives = None
        self.create_derivatives()
        self.der_encoder()
        self.der_regulation()
