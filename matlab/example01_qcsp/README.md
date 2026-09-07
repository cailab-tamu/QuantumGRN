# QuantumGRN in MATLAB — an annotated walkthrough

A MATLAB counterpart of `test/01_example.py`, written for reading rather than
for speed. Every quantum step goes through the **MATLAB Support Package for
Quantum Computing** (`quantumCircuit`, `ryGate`, `cryGate`, `simulate`,
`querystates`), and gradients come from the **parameter-shift rule** instead of
hand-derived derivative circuits.

It reproduces the Python result exactly: `max|θ_matlab − θ_python| = 1.1e-16`
after 100 epochs on the 6-gene × 28,208-cell LCL dataset.

## Run it

```matlab
cd matlab/example01_qcsp
check_conventions     % validate first — see below
qgrn_example01        % ~40 s
```

Requires the MATLAB Support Package for Quantum Computing (verified on R2026a
with support package 26.1.5).

## Which model is this?

This implements **Model A** — the circuit as published (Roman-Vicharra & Cai,
npj QI 9:67, 2023), i.e. the order-sensitive one. "Model A" / "Model B" is the
naming convention used in `test/order_free_rzz/` (defined in its
`Provenance.md` §0); Model B is the order-free RZZ variant proposed there, and
is **not** implemented here.

## The model in one paragraph

Each gene is one qubit. An encoder layer sets each qubit's marginal to that
gene's measured activation ratio, using `Ry(θ_kk)|0⟩ = cos(θ/2)|0⟩ + sin(θ/2)|1⟩`
so that `θ_kk = 2·asin(√activation)` reproduces it exactly. A regulation layer
then applies a controlled-`Ry(θ_kq)` from every gene *k* to every other gene *q*.
Measuring gives a distribution over the 2ⁿ activation patterns, which is fitted
by gradient descent to the pattern frequencies observed in the data. The fitted
off-diagonal angles **are** the network: magnitude is edge strength, sign is
activation vs repression.

## Files

| File | Role |
|---|---|
| `qgrn_example01.m` | main walkthrough — mirrors `01_example.py` step for step |
| `check_conventions.m` | validation; run this first |
| `qgrn_read_data.m` | binarize, order genes, build `p_obs` and activation ratios |
| `qgrn_circuit.m` | θ matrix → `quantumCircuit` |
| `qgrn_pout.m` | `simulate` → output distribution |
| `qgrn_laplace.m` | add-α smoothing (KL needs it; `p_obs` has exact zeros) |
| `qgrn_loss.m` | KL term + boundary constraint + squared error |
| `qgrn_dpraw.m` | parameter-shift derivatives |
| `qgrn_gradient.m` | full gradient via the chain rule |
| `qgrn_train.m` | gradient descent |
| `qgrn_draw_network.m` | network figure |
| `qgrn_pobs.m` | build `p_obs` for a given gene→qubit mapping |
| `ref_python_*.csv` | reference values exported from the Python package |

## Two things that will silently give you wrong answers

**1. Bit ordering.** The toolbox writes basis states so that character *k* is
qubit *k*: `ryGate(1,pi)` on 3 qubits gives `"100"`. The Python package indexes
states with gene 1 as the *least* significant bit. Getting this backwards
produces a perfectly plausible but completely wrong distribution — the mismatch
is 5.6e-01, not an obvious crash. `check_conventions.m` tests both orderings
against Python reference values; `q1lsb` is the correct one.

**2. The shift rule for controlled rotations.** The familiar two-term rule

```
df/dθ = [f(θ+π/2) − f(θ−π/2)] / 2
```

is exact only when the generator has eigenvalues ±1. For `c-Ry`, the generator
`|1⟩⟨1|⊗Y` has spectrum `{0, 0, +½, −½}`, so `f(θ)` carries **two** frequencies
and the two-term rule is wrong. The four-term rule is exact:

```
df/dθ = c₊[f(θ+π/2) − f(θ−π/2)] − c₋[f(θ+3π/2) − f(θ−3π/2)]
c± = (√2 ± 1) / (4√2)
```

The trap is that this is easy to miss: with a *single* controlled rotation on a
product state, `f` is a pure single-frequency sinusoid and the two-term rule
looks correct. Only when control branches interfere — i.e. in the real circuit —
does it break, and then by ~6%, which reads as slow convergence rather than as a
bug. `check_conventions.m` shows both cases.

Reference: Wierichs, Izaac, Wang & Lin, *General parameter-shift rules for
quantum gradients*, Quantum **6**, 677 (2022).

## This model is order-sensitive

Genes are mapped to qubits **sorted by activation ratio** (`qgrn_read_data.m`),
and regulation layers are emitted one control at a time (`qgrn_circuit.m`).
`c-Ry` layers with different controls do not commute: applying `c-Ry(k→p)` then
`c-Ry(p→j)` is not the same circuit as the reverse, because the second gate's
control has already been rotated. So the fitted θ depends on the gene ordering,
while the symmetrised gradient imposes `θ_kp = θ_pk` on a structurally
asymmetric circuit. This is the published criticism recorded in
`test/order_free_rzz/Provenance.md`.

It has been measured, using the model code in this folder. The scripts live
with the ordering analysis they feed, in `test/order_free_rzz/`:
`ordering_sweep.m` (5 orderings, 100 epochs, MATLAB, ~3 min) and
`big_sweep.py` (41 orderings × 1200 epochs, Python for speed, ~20 min).
Recorded as finding **F6** in that folder's `Provenance.md`.

```
paper ordering KL = 4.250449e-03
random  KL: min 3.650e-03  median 4.105e-03  max 4.467e-03
paper z-score vs random: +0.61      percentile 62.5%

max edge spread 0.0162 (IRF4-PRDM1)   max std 0.0039   sign flips 0 of 15
```

Three things follow.

**The effect is real but second-order.** Every fitted `|θ|` is well below π/4, so
the `c-Ry` gates nearly commute. This confirms and tightens finding F2 of
`Provenance.md` (spread 0.015, std 0.005, all signs preserved, 6 orderings) at
41 orderings and full convergence. Relative spread is 8–19% for the strong
edges.

**The paper's ordering is unremarkable.** It sits at the 62.5th percentile of
the random distribution (z = +0.61); 15 of 40 random orderings fit worse. A
5-ordering pilot suggested the paper ordering was the *worst* of its set — that
did not survive the larger sweep and was an artifact of four random draws.

**But edge *presence* is not stable at the pruning boundary.** The readout
prunes at `|θ| > π/180 ≈ 0.01745` (`qgrn_draw_network.m`, matching
`qscgrn/visualization/qsc_grn.py:46`). Across 41 orderings nine edges clear it
every time and four never do, but two straddle it:

| edge | paper θ | present in |
|---|---|---|
| IRF4–AICDA | −0.0182 | 24 / 41 |
| REL–RELA | +0.0147 | 4 / 41 |

So "all signs preserved" holds while the *edge set* does not. IRF4–AICDA is the
edge the paper flags as incorrectly recovered — whether it appears at all is a
function of gene ordering.

If you need order-independence rather than a bound on the error, the fix is not
in this folder — `test/order_free_rzz` replaces the `c-Ry` cascade with a
commuting RZZ layer (Model B), which has no ordering at all and fits better
(KL 1.74e-3).

## One inherited approximation, kept on purpose

The published model builds `dp/dθ` from the **raw** Born probabilities while
evaluating the loss on the `|0…0⟩`-rescaled ones — it drops the derivative of
the rescaling step. This port reproduces that choice so the two implementations
agree numerically, and `check_conventions.m` §4 quantifies the gap: on a test
case where `|0…0⟩` holds 57% of the raw mass, the model gradient is `+0.636`
where the true derivative of the stated loss is `+1.465`.

This is worth knowing if you extend the model. It is not a bug to fix silently —
changing it changes the published method.

## Cost

Parameter-shift needs 4 circuit simulations per controlled angle per epoch:
`n²−n = 30` angles × 4 × 100 epochs ≈ 12,000 simulations of a 6-qubit,
36-gate circuit, about 40 s. The Python package's hand-derived derivative
circuits are cheaper; parameter-shift is used here because it generalises to any
gate set and is the technique that carries over to hardware.
