# Isomorphic Engine Experiments

Seven Rust experiments verify the Yang-Mills mass gap using the Isomorphic Engine (12 parallel CPU solvers, sparse CSR pipeline, tested to N = 8,232 spins).

## Experiments

| File | Purpose | Runtime |
|------|---------|---------|
| `yang_mills_mass_gap.rs` | Killing form identity, U_24 membership, basic lattice | 15s |
| `yang_mills_confinement.rs` | Center symmetry, barrier scaling B(L), topology | 94s |
| `yang_mills_chaos.rs` | Lyapunov exponents, GUE KS convergence | 4s |
| `yang_mills_continuum.rs` | Multi-beta mass gap, asymptotic freedom | 3s |
| `yang_mills_high_dim.rs` | Scale survey to N=5,184, barrier to N=8,232 | 3.3hr |
| `yang_mills_os_axioms.rs` | f(L) convergence, C(r) decay, OS axioms | 147s |
| `yang_mills_falsifiability.rs` | 15/15 falsifiable predictions | 183s |

## Running

Copy experiments to an Isomorphic Engine installation:

```bash
cp *.rs /path/to/isomorphic-engine/examples/
cd /path/to/isomorphic-engine
cargo run --release --example yang_mills_falsifiability
```

## Key Results

- **Barrier**: B(L) ~ L^3.18 (SU(2)), L^3.09 (SU(3))
- **Chaos**: lambda_max = 0.19-0.28 (9 configurations, all positive)
- **GUE**: KS distance 0.68 -> 0.22 (decreasing with N)
- **Mass gap**: Delta > 0 at all 24 tested configurations
- **OS axioms**: f(L) converges < 1%, C(r) ~ exp(-2.3r)
