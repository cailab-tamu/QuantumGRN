# NEXT_STEPS.md — open items (owner: James; last updated 2026-09-06)

Priority order. Tick and date when done.

- [ ] **Run `matlab/check_qcsp.m` in MATLAB.** Confirms plain code vs
      Quantum Computing Support Package agree and fixes the basis-string bit
      order if needed. Produces `results/ising_bm_mix_pi6.qasm`.
- [ ] **Ordering sweep at n = 10–15 genes** (Model A) to find where the
      reviewer's criticism becomes material. Needs a gene-selection rule
      (highly variable genes) and the merged LCL matrix from scGEAToolbox.
- [ ] **Mixer angle**: make θ_mix trainable, or sweep finely in [π/8, π/3];
      π/3 already degrades (KL 8.6e-3, one sign flip).
- [ ] **Identifiability of θ in Model B**: at minimum characterise the
      RZZ(θ+π) aliasing under the mixer; keep L_cons. Possibly derive the
      closed-form p_out (Walsh–Hadamard of phased amplitudes at mixer π/2).
- [ ] **Redo the boundary-effect analysis** (paper Fig. 4) for RZZ + Rx.
- [ ] **Directed variant** (option #3 from the discussion): drop
      symmetrisation in Model A and treat the best-fitting ordering as an
      inferred regulatory hierarchy. Independent of Model B.
- [ ] **Analytic gradients** (parameter-shift for RZZ/Rx) to replace finite
      differences in the MATLAB port before scaling n.
- [ ] **Integrate Model B into `qscgrn`** as an alternative circuit type once
      the above are settled.
- [ ] **Write-up**: rebuttal paragraph (finding F2) + new-method section
      (findings F3–F5) for a response letter or follow-up paper.

## Done
- 2026-09-05  Python prototype, both models, reproduction of Fig. 3c.
- 2026-09-06  MATLAB/Octave port, verified; support-package variants written.
