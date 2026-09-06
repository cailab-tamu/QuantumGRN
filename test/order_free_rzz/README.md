# order_free_rzz — order-free qscGRN (Ising Born machine)

Exploratory prototype addressing the gene-ordering sensitivity of qscGRN
(Roman-Vicharra & Cai, npj QI 2023) by replacing the sequential c-Ry
regulation layers with a commuting all-pairs RZZ layer plus an Rx mixer.

**Start here:** `Provenance.md` (what, why, findings, file map).
**For Claude Code sessions:** `CLAUDE.md`. **Open items:** `NEXT_STEPS.md`.

Quick run (Python, ~20 min):
```
pip install autograd numpy pandas
cd test/order_free_rzz && python run_comparison.py
```
Quick run (MATLAB/Octave, ~30 min in Octave):
```
cd test/order_free_rzz/matlab && matlab -batch run_comparison
```
Headline: Model A reproduces Fig. 3c exactly (KL 4.25e-3); ordering changes
theta by <= 0.015 with all signs preserved; Model B (RZZ + Rx(pi/6), diagonals
trained) reaches KL 1.74e-3 with 10/10 sign agreement and no ordering at all.
