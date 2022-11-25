import numpy as np
from ..utils import info_print


__all__ = ["qscgrn_model"]


class qscgrn_model():
    """
    Construct and executes a quantum circuit for the QuantumGRN model.
    Attributes
    ----------
    ngenes : int
        Number of genes in the Quantum GRN.
    genes : list
        Gene list for QuantumGRN modelling.
    theta : pd.Series
        Theta values given edges in the Quantum GRN.
    edges : list
        Edges for the QuantumGRN.
    drop_zero : bool
        If True, a normalization is applied to the probability distribution
        in the state `|0>_n`, and the rest of the distribution is rescale.
    Methods
    -------
    _enc_layer()
        Constructs the encoder layer of the quantum circuit.
    _reg_layer()
        Constructs the regulation layer of the quantum circuit.
    _meas_layer()
        Construct the measument layer of the quantum circuit.
    _compile_run()
        Compiles and runs the simulation of the quantum circuit.
    _qiskit_2_np()
        Process the output of the quantum circuit to return a np.array.
    run_qiskit()
        Executes methods to return the probability distribution.
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
        self.ngenes =len(genes)
        self.genes = genes
        self.theta = theta
        self.edges = edges
        self.drop_zero = drop_zero


    def __str__(self):
        return ("Quantum circuit for {ngenes} using Qiskit package"
                "for GRN modelling".format(ngenes=len(self.genes)))


    def _enc_layer(self, circ):
        """
        Constructs the encoder layer. Reset the input register to `0` and
        initialize the `R_y` gates to the corresponding values.
        Parameters
        ----------
        circ : qiskit.QuantumCircuit
            The quantum circuit to append the encoder layer `L_enc`.
        Returns
        -------
        circ : qiskit.QuantumCircuit
            A quantum circuit with the `L_enc` layer appended.
        """
        for i, gene in enumerate(self.genes):
            circ.reset(i)
            circ.ry(self.theta[(gene, gene)], i)

        circ.barrier(range(self.ngenes))
        return circ


    def _reg_layer(self, circ, threshold):
        """
        Constructs the regulation layers with the corresponding
        values for each controlled `c-R_y` gate.
        Parameters
        ----------
        circ : qiskit.QuantumCircuit
            The quantum circuit to append the encoder layer `L_enc`.
        threshold : float
            The threshold for dropping edges from the network and should be
            the same than the qscgrn.visualization.qsc_grn.draw_network method.
            Details in the manuscript.
        Returns
        -------
        circ : qiskit.QuantumCircuit
            A quantum circuit with the `L_reg` layers appended.
        """
        from ..qcircuit import edges_to_index

        theta = self.theta[self.edges]
        theta = theta[np.abs(theta) > (threshold * np.pi / 180)]
        if len(theta) == 0:
            return circ

        edges = theta.index.to_list()
        idxs = edges_to_index(self.genes, edges)

        flag = idxs[0][0]
        for idx, edge in zip(idxs, edges):
            if idx[0] != flag:
                circ.barrier(range(self.ngenes))
                flag = idx[0]
            circ.cry(theta[edge], idx[0], idx[1])

        return circ


    def _meas_layer(self):
        """
        Constructs the measurement layer of the quantum circuit.
        Returns
        -------
        meas : qiskit.QuantumCircuit
        """
        from qiskit import QuantumCircuit
        meas = QuantumCircuit(self.ngenes, self.ngenes)
        meas.barrier(range(self.ngenes))
        meas.measure(range(self.ngenes), range(self.ngenes))
        return meas


    def _compile_run(self, qc):
        """
        Compiles and runs the quantum circuit. Currently using only
        the Aer simulator.
        Parameters
        ----------
        qc : qiskit.QuantumCircuit
            The quantum circuit for the GRN model.
        Returns
        -------
        counts : qiskit.counts
            The counts distribution for each state in Qiskit format.
        """
        from qiskit import transpile
        # the next line can be changed to use different simulators
        from qiskit.providers.aer import AerSimulator

        backend = AerSimulator()
        qc_compiled = transpile(qc, backend)
        job_sim = backend.run(qc_compiled, shots=8192)
        result_sim = job_sim.result()
        counts = result_sim.get_counts(qc_compiled)
        return counts


    def _qiskit_2_np(self, counts, drop_zero):
        """
        Translates the Qiskit format to a readable array.
        Parameters
        ----------
        counts : qiskit.counts
            The count distribution in Qiskit format.
        drop_zero : bool
            If True, a normalization step is computed.
        Returns
        -------
        prob : np.array
            The count distribution in a readable format.
        """
        prob = np.zeros((2**self.ngenes,))
        for idx in counts:
            prob[int(idx, 2)] = counts[idx]

        if drop_zero:
            prob[0] = 0.
            return prob / np.sum(prob)
        else:
            return prob


    def run_qiskit(self, threshold=1, filename=None):
        """
        Joins the layers and executes the quantum circuit
        Parameters
        ----------
        threshold : float
            The threshold for dropping edges from the network and should be
            the same than the qscgrn.visualization.qsc_grn.draw_network method.
            Details in the manuscript.
        filename : str
            The file name where the schematic of the quantum circuit is saved,
            given a filename.
        Returns
        -------
        counts : np.array
            The count distribution in a readable format.
        """
        from qiskit import QuantumCircuit
        circ = QuantumCircuit(self.ngenes)
        circ = self._enc_layer(circ)
        circ = self._reg_layer(circ, threshold)
        meas = self._meas_layer()
        qc = meas.compose(circ, range(self.ngenes), front=True)
        counts = self._compile_run(qc)

        if filename is not None:
            qc.draw(output="mpl", filename=filename)
            info_print("Drawing the quantum circuit of the qscGRN model "
                       "and saving to {file}".format(file=filename))

        return self._qiskit_2_np(counts, self.drop_zero)

