# Contributing and Reproducibility Guide

## Quick Start

### Python Verification (59/59 checks, ~5 minutes)

```bash
# Option A: Conda
conda create -n ym-verify python=3.10 numpy scipy matplotlib
conda activate ym-verify
python scripts/verify_yang_mills.py

# Option B: pip
python -m venv .venv && source .venv/bin/activate  # or .venv\Scripts\activate on Windows
pip install numpy scipy matplotlib
python scripts/verify_yang_mills.py
```

### Figure Regeneration

```bash
python scripts/generate_figures.py   # Produces 8 PNGs in figures/
```

### Engine Experiments (Rust)

The engine experiments require the Isomorphic Engine (not included). The Rust source code for all 7 experiments is in `engine/`. To run with an existing engine installation:

```bash
# Copy experiment files to your engine's examples/ directory
cp engine/*.rs /path/to/isomorphic-engine/examples/

# Build and run
cd /path/to/isomorphic-engine
cargo run --release --example yang_mills_mass_gap
cargo run --release --example yang_mills_confinement
cargo run --release --example yang_mills_chaos
cargo run --release --example yang_mills_falsifiability
```

## What You Can Verify Independently

- **Killing form identity** (exact algebra): Compute f^{abc} for SU(3), form J^{ab} = sum f^{acd} f^{bcd}, check J = 3 I_8 and Tr = 24.
- **SU(3) uniqueness**: Check N(N^2-1) = 24 has no solution for N != 3.
- **E_8 root system**: Construct 240 roots, verify Cartan matrix and subgroup Casimirs.
- **GUE statistics**: Generate GUE random matrices, compute KS distance vs Wigner surmise.
- **Asymptotic freedom**: Verify beta-function formula with C_2(SU(3)) = 3.

## For Reviewers

1. **Start with** `scripts/verify_yang_mills.py` — runs all 59 checks in ~5 minutes.
2. **Read** `PROOF.md` for the complete self-contained proof.
3. **Check** `data/verification-summary/verification_summary.json` for detailed check-by-check results.
4. **Examine** `papers/yang-mills/main.tex` for the full paper with all theorems and proofs.

## Reporting Issues

Please open a GitHub issue for any:
- Failed verification checks
- Mathematical errors in the proof
- Computational irreproducibility
- Missing or incorrect data
