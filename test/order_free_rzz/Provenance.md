# Provenance — order-free qscGRN via commuting RZZ layer (Ising Born machine)

**Date:** 2026-09-05
**Author of analysis:** Claude (Anthropic), in conversation with James J. Cai
**Repo context:** `cailab-tamu/QuantumGRN`, branch `main`, shallow clone taken 2026-09-05
**Status:** exploratory prototype; not integrated into the `qscgrn` package

---

## 1. Motivation

Source paper: Roman-Vicharra C, Cai JJ. *Quantum gene regulatory networks.*
npj Quantum Information 9:67 (2023). https://doi.org/10.1038/s41534-023-00740-6

A peer paper (citing qscGRN as its ref. [66]) raised the criticism:

> "the circuit construction is sensitive to gene ordering, as it prioritizes
> genes with higher activation ratios in the qubit mapping."

Diagnosis (see §3): the criticism is correct in principle. Within a regulation
layer L_k all c-R_y gates share control k and commute; across layers,
c-R_y(k→p) and c-R_y(p→j) do not, because the second gate's control has already
been rotated. Eq. 4 of the paper is therefore a *sequential* model in which
g_0 (highest activation ratio) acts first. The symmetrised gradient (Eq. 8)
imposes θ_kp = θ_pk on a circuit that is structurally asymmetric.

Four remedies were discussed; James selected **#4, an order-free ansatz**:

1. Empirical sensitivity bound (refit under random permutations).
2. Permutation marginalisation during training.
3. Drop symmetrisation; treat the ordering as a learnable regulatory hierarchy.
4. Replace the gate sequence with a symmetric generator so no ordering exists.

## 2. Source of the design for #4

Miller L, Lee Y. *QuIC: A Training-Free Quantum Graph Embedding from Ideal
Analysis to Practical Hardware Evaluation.* arXiv:2604.18841v1 (20 Apr 2026).

Only the **entangling layer** was borrowed: a product of RZZ gates, one per
edge. RZZ gates are diagonal in the computational basis, hence commute, hence
the layer is order-free and symmetric in (k,p). With trainable angles,
U_ent(θ) = exp(−i Σ_{k<p} θ_kp Z_k Z_p / 2) **exactly** (no Trotter error).

Not borrowed from QuIC, and why:

| QuIC component | Used? | Reason |
|---|---|---|
| RZZ entangler, commuting | yes | the order-free coupling layer |
| RX mixer after entangler | yes (angle re-tuned) | diagonal phases are invisible in Z-basis without it |
| Degree-proportional encoder | no | replaced by activation-ratio R_y as in qscGRN |
| Fixed / training-free parameters | no | θ must be fitted to p_obs |
| Sorted output distribution | no | destroys gene labels |
| Injectivity theorem (Thm 3) | no | concerns statevectors of unweighted graphs, not identifiability of continuous θ from measured p_out (their Remark 10) |

The resulting model is an **IQP-type Ising Born machine**; cf. Coyle B et al.,
"The Born supremacy: quantum advantage and training of an Ising Born machine",
npj Quantum Information 6:60 (2020).

## 3. Models implemented

Both in `qscgrn_vs_ising_bm.py`, exact statevector simulation, n = 6 qubits,
q0 = most significant bit of the label.

**Model A — original qscGRN (reference).**
|0⟩^n → R_y(θ_kk) on each qubit → for k = 0..n−1, for p ≠ k: c-R_y(θ_kp) with
control k, target p → measure. Diagonals frozen at 2·asin(√act_k).
Faithful to Algorithms 1–3 of the paper.

**Model B — Ising Born machine (proposed).**
|0⟩^n → R_y(θ_kk) on each qubit → Π_{k<p} RZZ(θ_kp) → R_x(θ_mix) on each qubit
→ measure. Diagonals **trained** (see finding F3). θ_mix ∈ {π/6, π/4, π/3}
swept; not trained.

**Common.** p_obs and p_out with |0⟩^n set to 0 and rescaled; Laplace
smoothing α = 1; loss = KL(p̂_out ‖ p̂_obs) + λ·L_cons with the paper's dynamic λ;
plain gradient descent, lr = 0.05, gradient symmetrised (G + Gᵀ)/2; 1200
iterations, no early stop. Gradients by reverse-mode autodiff (`autograd`).

## 4. Data

`dataset/expr_matrix_pearsonresidual_7.txt` from the repo: Pearson-residual
normalised expression of IRF4, REL, PAX5, RELA, PRDM1, AICDA in 28,208 LCL
cells (GSE126321 + GSE158275). Binarised at 0. Activation ratios: IRF4 0.349,
REL 0.345, PAX5 0.330, RELA 0.311, PRDM1 0.272, AICDA 0.207 (this is the paper
ordering).

## 5. Findings (full numbers in `results/results.log`)

**F1 — Reproduction.** Model A at 1200 iterations gives KL = 4.25e−3 and
reproduces the Fig. 3c adjacency matrix to three decimals
(IRF4–PRDM1 0.166, REL–PRDM1 −0.129, PAX5–PRDM1 −0.101, REL–AICDA 0.084,
REL–PAX5 0.078, PAX5–AICDA 0.067, PAX5–RELA 0.043, IRF4–REL −0.070,
IRF4–AICDA −0.018, PRDM1–AICDA −0.078).

**F2 — Ordering sensitivity is real but small in this regime.** Six random
gene orderings, mapped back to a common index: maximum spread of any θ_kp is
0.015 (IRF4–PRDM1), maximum std 0.005, all 15 signs preserved. Because
fitted |θ| < π/4, the c-R_y gates nearly commute and the effect is
second-order. This is a quantitative rebuttal to the reviewer's paragraph for
the published 6-gene analysis; it says nothing about larger n or larger θ.

**F3 — Frozen diagonals fail in Model B.** With θ_kk fixed, KL stalls at
~3e−2. After an entangler plus mixer the single-qubit marginals depend on the
couplings, so the encoder no longer pins act_k. Diagonals must be trained;
they drift (e.g. PRDM1 1.098 → 0.939, AICDA 0.945 → 0.841 at θ_mix = π/6).

**F4 — Model B works and fits better.** θ_mix = π/6, diagonals trained:
KL = 1.74e−3 (vs 4.25e−3 for A); 10/10 sign agreement with Fig. 3c;
edge magnitudes ~1.3–1.5× larger. IRF4–AICDA — the edge the paper reports as
incorrectly recovered — comes out at −0.003, below the pruning threshold
π/(180·2), i.e. absent rather than wrong-signed. No gene ordering exists in
the model.

**F5 — Mixer angle.** π/6 and π/4 behave similarly (KL 1.7e−3, 2.9e−3);
π/3 degrades (KL 8.6e−3, PRDM1–AICDA flips sign, 9/10 concordance).
Suggested default π/6–π/4; could be made trainable.

## 6. Open issues / next steps

- Identifiability of θ from p_out in Model B is not established. Keep the
  L_cons constraint: RZZ(θ+π) differs from RZZ(θ) only by Z⊗Z, which aliases
  under the mixer.
- The boundary-effect analysis of the paper (Fig. 4) is specific to c-R_y and
  would need redoing for RZZ + R_x.
- Repeat F2 at n = 10–15 to find where ordering sensitivity actually bites.
- Compare against the earlier tensor-network / Born-machine work in the lab;
  the RZZ layer is the commuting special case of that family.
- Hardware: all-to-all RZZ is n(n−1)/2 two-qubit gates; trivially within depth
  budgets at n ≤ 15 (QuIC reports a ~210–250-layer limit on ibm_fez).

## 7. Files

```
test/order_free_rzz/
├── Provenance.md                       this file
├── README.md                           one-screen orientation
├── CLAUDE.md                           conventions + status for AI coding sessions
├── NEXT_STEPS.md                       open items and done log
├── qscgrn_vs_ising_bm.py               models, loss, fit(); run as script for early-stop sweep
├── run_comparison.py                   fixed-1200-iteration comparison; source of results.log
├── matlab/                             MATLAB/Octave port (plain statevector code, no toolbox needed)
│   ├── run_comparison.m                driver; same experiments as run_comparison.py
│   ├── pout_qscgrn.m                   Model A forward model
│   ├── pout_ising_bm.m                 Model B forward model
│   ├── qgrn_fit.m                      GD with symmetrised finite-difference gradient
│   ├── qgrn_loss.m, build_pobs.m, read_expr.m, apply_1q.m, apply_c1q.m,
│   │   gate_defs.m, zz_signs.m, print_theta.m
│   ├── pout_qscgrn_qcsp.m, pout_ising_bm_qcsp.m, qcsp_probs.m
│   │                                   same models built with the Quantum Computing Support Package
│   └── check_qcsp.m                    cross-check plain code vs support package; exports QASM (MATLAB only, UNTESTED)
└── results/
    ├── octave_run.log                  console output of run_comparison.m (Octave 8.4)
    ├── results_matlab.mat              thA, stack (4 orderings), thB{1..3}, mixes, names, paper
    ├── results.log                     console output of run_comparison.py
    ├── results_earlystop.log           output of qscgrn_vs_ising_bm.py (early-stop run, 8 orderings)
    ├── thetaA_paper_order_1200it.npy   Model A, paper ordering (= Fig. 3c)
    ├── thetaA_8_orderings_earlystop.npy  Model A, 8 random orderings, early-stop run (index order = paper order)
    ├── thetaB_mix_pi6.npy              Model B, θ_mix = π/6, diagonals trained
    ├── thetaB_mix_pi4.npy              Model B, θ_mix = π/4
    └── thetaB_mix_pi3.npy              Model B, θ_mix = π/3
```

Run from `test/order_free_rzz/` (data path is `../../dataset/…`):

```
pip install autograd numpy pandas
python run_comparison.py > results/results.log
```

## 7b. MATLAB port

`matlab/run_comparison.m` reproduces the Python experiments with no dependence
on the Quantum Computing Support Package (unlike `matlab/i_fullcirc_*.m` in
the repo root). Gradients are central finite differences (h = 1e-6) on the
symmetric parameterisation, which equals the paper's (∇+∇ᵀ)/2 update.
Tested in Octave 8.4; MATLAB-compatible syntax throughout. Results agree with
the Python run to three decimals: Model A KL 4.25e-3 and Fig. 3c reproduced;
max θ spread across 4 orderings 0.013; Model B π/6 KL 1.74e-3, 10/10 signs;
π/4 2.89e-3, 10/10; π/3 8.61e-3, 9/10. Runtime ≈ 4 min per fit in Octave
(finite differences; ~8 fits total); MATLAB will be faster. Run from
`test/order_free_rzz/matlab/`.

**Support-package cross-check (untested).** `check_qcsp.m` rebuilds both
models with `ryGate`/`cryGate`/`rzzGate`/`rxGate`, simulates with the
support package, and asserts agreement with the plain code at the fitted θ
(tries both basis-string bit orderings and reports which matches). It also
writes `results/ising_bm_mix_pi6.qasm` via `generateQASM` for hardware runs.
This part requires MATLAB and could not be run in the analysis environment;
the plain-code results above do not depend on it.

## 8. Environment

Python 3.12, numpy 2.4.4, pandas 3.0.2, autograd (latest PyPI, 2026-09).
Exact statevector simulation only; no Qiskit, no shot noise, no hardware.
Runtime ≈ 20 min for `run_comparison.py` on one CPU core.
