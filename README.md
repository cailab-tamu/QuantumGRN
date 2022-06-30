# QuantumGRN
A quantum circuit model for inferring gene regulatory networks with application to single-cell transcriptomic data

Quantum computing holds the promise to achieve certain types of computation that would otherwise be unachievable by classical computers. The advent in the development of quantum algorithms has enabled a variety of applications in chemistry, finance, and cryptography. In this work, we present a quantum circuit model for constructing gene regulatory networks (GRNs) from single-cell transcriptomic data. The model is based on the idea of using qubit-qubit entanglement to simulate the interactions between genes. Each qubit in the circuit represents a gene, and qubits are entangled to simulate the interaction between genes. The strength of gene interactions is estimated using the rotation angle of controlled unitary gates between qubits. We provide preliminary results that suggest our quantum single-cell GRN (qscGRN) modeling method is competitive and warrants further investigation. Specifically, we present the preliminary results derived from the single-cell RNA sequencing (scRNA-seq) data of human lymphoblastoid cell lines, focusing on genes in the nuclear factor-kappa B (NF-κB) signaling pathway. We demonstrate that our qscGRN model can recover known and detect novel regulatory relationships, setting the stage for further investigations on GRNs, given that relationships between fully interconnected genes are approached more effectively by quantum modeling than by statistical correlations. Our quantum circuit model enables the modeling of vast feature space occupied by cells in different transcriptionally activating states, simultaneously tracking activities of thousands of interacting genes and constructing more realistic single-cell GRNs without relying on statistical correlation or regression. We anticipate that quantum computing algorithms based on our circuit model will find more applications in data-driven life sciences, paving the way for the future of predictive biology and precision medicine.

## Installation (Coming soon...)
Use the package manage [pip](https://pip.pypa.io/en/stable/) to install QuantumGRN
```bash
pip install QuantumGRN
```

## Run example
Run the example.py file using QuantumGRN package
```bash
cd test
python example.py
```

## Contributing
Pull requests are welcome.

## License
[MIT](https://choosealicense.com/licenses/mit/)