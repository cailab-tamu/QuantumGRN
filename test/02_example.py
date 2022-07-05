import os
import numpy as np
import pandas as pd
from qscgrn import *

path = os.getcwd()
filename="../dataset/expr_matrix_pearsonresidual_7.txt"
df = pd.read_csv(os.path.join(path, filename), delimiter='\t')
df = df.set_index('genes').T

ncells, ngenes = df.shape
df = qsc_order_gene(df)
genes = df.columns.to_list()
p_obs = qsc_distribution(df)
activation = qsc_activation_ratios(df)

mask = mini_hist(ngenes, p_obs, limit=0.01, ymax=0.15,
                 title="Observed distribution",
                 filename="02_p_obs.svg")

theta = theta_init(genes, activation_ratios=activation)
edges = edges_init(genes)
qgrn = model(ncells, genes, theta, edges, p_obs, epochs=1000)
qgrn.train()
p_out =qgrn.p_out.reshape(2**ngenes,)
p_qiskit = qgrn.run_qiskit(filename="02_qgrn_circuit.svg")

comparison_hist(ngenes, p_obs, p_out, limit=0.01, ymax=0.15, mask=mask,
                filename="02_comparison_obs_out.svg", tags=["p^{obs}", "p^{out}"])

comparison_hist(ngenes, p_obs, p_qiskit, limit=0.01, ymax=0.15, mask=mask,
                filename="02_comparison_obs_qiskit.svg", tags=["p^{obs}", "p^{qiskit}"])

draw_network(genes, edges, qgrn.theta, filename="02_qgrn_network.svg")
