# Data Dictionary

All data files are JSON format. No external downloads required.

## Coupling Matrices

| File | Records | Description |
|------|---------|-------------|
| `coupling-matrices/coupling_matrix_JG.json` | SU(2), SU(3), SU(2)-SU(8) | Yang-Mills coupling matrices J_G, eigenvalues, traces |

Key fields: `J` (matrix), `eigenvalues`, `trace`, `C2_adj`, `dim_g`, `trace_equals_omega`

## E_8 Root System

| File | Records | Description |
|------|---------|-------------|
| `e8-root-system/e8_root_system.json` | 240 roots | Root count, minimal norm, Casimir traces for all subgroups |

## Barrier Scaling

| File | Records | Description |
|------|---------|-------------|
| `barrier-scaling/lattice_results.json` | SU(2), SU(3) | Lattice gauge theory: plaquettes, Wilson loops, mass gap estimates |

Engine results (from `yang_mills_high_dim.rs`):
- SU(2): 6 data points, L = 3-8, N up to 4,608, alpha = 3.18
- SU(3): 5 data points, L = 3-7, N up to 8,232, alpha = 3.09

## Mass Gap

| File | Records | Description |
|------|---------|-------------|
| `mass-gap/mass_gap_bounds.json` | 3 groups | Explicit bounds Delta >= C_2 Lambda / sqrt(24) |

Engine results: Delta > 0 at all 24 tested (G, L, beta) configurations.

## GUE Statistics

| File | Records | Description |
|------|---------|-------------|
| `gue-statistics/spectral_statistics.json` | 200 eigenvalues | KS test, spacing variance, GUE metrics |

## OS Axioms

Engine results from `yang_mills_os_axioms.rs`:
- f(L) = E_0/V converges to < 1% (SU(2), L = 7-8)
- C(r) ~ exp(-2.30 r), Delta_eff = 2.30
- |C(1, L=6) - C(1, L=7)| = 0.004

## Leech Lattice

| File | Records | Description |
|------|---------|-------------|
| `leech-lattice/leech_lattice.json` | 1 | dim = 24, min_norm = 2, kissing = 196,560 |

## Verification Summary

| File | Records | Description |
|------|---------|-------------|
| `verification-summary/verification_summary.json` | 59 checks | All check results with pass/fail and details |

Central result: `Tr(J_YM^{SU(3)}) = 24 = Omega`
