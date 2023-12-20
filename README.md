# QuantumGRN
A quantum circuit model for inferring gene regulatory networks with application to single-cell transcriptomic data. Reference in the [manuscript](https://www.nature.com/articles/s41534-023-00740-6).

we introduce a quantum circuit model for inferring gene regulatory networks (GRNs) from single-cell transcriptomic data. The model employs qubit entanglement to simulate interactions between genes, resulting in competitive performance and promising potential for further exploration. We applied our quantum GRN modeling approach to single-cell transcriptomic data from human lymphoblastoid cells, focusing on a small set of genes involved in innate immunity regulation. Our quantum circuit model successfully predicted the presence and absence of regulatory interactions between genes, while also estimating the strength of these interactions. We argue that the application of quantum computing in biology has the potential to provide a better understanding of single-cell GRNs by more effectively approaching the relationship between fully interconnected genes compared to conventional statistical methods such as correlation and regression. Our results encourage further investigation into the creation of quantum algorithms that utilize single-cell data, paving the way for future research into the intersection of quantum computing and biology.

## Installation (test Pypi)
Use the package manage [pip](https://pip.pypa.io/en/stable/) to install QuantumGRN
```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple QuantumGRN
```

## Run example
Run the examples file using QuantumGRN package
```bash
cd test
python 02_example.py
```

## Usage
QuantumGRN is used as normal package,
```python
import numpy as numpy
import pandas as pd
from qscgrn import *

df = pd.read_csv("dataset/expr_matrix_pearsonresidual_7.txt", delimiter='\t')
df = df.set_index('genes').T

ncells, ngenes = df.shape
df = qsc_order_gene(df)
genes = df.columns.to_list()
p_obs = qsc_distribution(df)
activation = qsc_activation_ratios(df)

theta = theta_init(genes, activation_ratios=activation)
edges = edges_init(genes)
qgrn = model(ncells, genes, theta, edges, p_obs)
qgrn.train()

draw_network(genes, edges, qgrn.theta, filename="qgrn_network.png")
```

More detailed code in the tutorials directory

## Contributing
Pull requests are welcome.

## License
[MIT](https://choosealicense.com/licenses/mit/)
