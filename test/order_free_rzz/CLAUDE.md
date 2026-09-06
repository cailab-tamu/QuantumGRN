# CLAUDE.md — context for AI coding sessions in test/order_free_rzz

Read `Provenance.md` first. This file records conventions and status so a
fresh session can continue without the original chat.

## What this folder is
A self-contained prototype comparing
- **Model A**: the published qscGRN circuit (Ry encoder → sequential c-Ry
  layers). Order-sensitive by construction.
- **Model B**: Ry encoder → all-pairs RZZ(θ_kp) → Rx(mixer) → measure.
  Order-free (RZZ gates commute). An IQP-type Ising Born machine.
Both fitted to the 6-gene LCL p_obs with the paper's KL + L_cons loss.
Python (`autograd`) and MATLAB/Octave (finite differences) implementations
give identical results to 3 decimals.

## Conventions — do not change silently
- Gates: Ry(t) = exp(-i t/2 Y), Rx(t) = exp(-i t/2 X), RZZ(t) = exp(-i t/2 Z⊗Z)
  (Qiskit conventions). Model B applies the RZZ phase as a diagonal vector.
- State indexing: **Python** q0 = most significant bit; **MATLAB** qubit 1 =
  least significant bit (column-major). θ matrices are identical either way.
- |0…0⟩ probability is zeroed and the rest rescaled in both p_obs and p_out.
- Laplace smoothing α = 1 with m = 28,208 cells; dynamic λ = L_KL / L_cons
  recomputed every iteration; lr = 0.05; gradient symmetrised (G+Gᵀ)/2;
  1200 iterations, no early stop (the paper's threshold 2ⁿ·1e-4 stops A at
  ~60 iterations with visibly smaller θ — see results/results_earlystop.log).
- Gene order = descending activation ratio: IRF4 REL PAX5 RELA PRDM1 AICDA.
- Model B diagonals **must be trained**; frozen diagonals stall at KL ~3e-2.
- Pruning threshold from the paper: |θ| < π/(180·2) ≈ 0.0087 → no edge.

## Status
| Item | State |
|---|---|
| Python A/B fits, ordering sweep | done, verified (results/results.log) |
| MATLAB port (plain statevector) | done, verified in Octave 8.4 (results/octave_run.log) |
| Support-package variants + `check_qcsp.m` + QASM export | written, **untested** (needs MATLAB + Quantum Computing Support Package) |
| Integration into `qscgrn` package | not started |
| n > 6 experiments | not started |

## Data
`../../dataset/expr_matrix_pearsonresidual_7.txt` (repo file; 6 × 28,208).
Scripts use relative paths; run from this folder (Python) or `matlab/` (MATLAB).

## Working rules
- Keep `Provenance.md` §5 (findings) and this Status table current when
  results change. Add a dated line to `NEXT_STEPS.md` when an item closes.
- Do not modify the `qscgrn` package or `matlab/i_fullcirc_*.m` in the repo
  root from here; propose changes first.
- When adding a new circuit variant, add it as a new `pout_*.py/.m` function
  and re-run the Fig. 3c sign-concordance check (10/10 for A and B at π/6).
- Reference for the borrowed entangler: Miller & Lee, QuIC, arXiv:2604.18841.
  Reference for the model class: Coyle et al., npj QI 6:60 (2020).
