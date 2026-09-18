# HCN···HCl, dense chlorine grid variant (`cl400`)

An alternate copy of the HCN···HCl set in `../hcn_hcl/` that differs in one
setting: 400 radial shells on chlorine instead of the 110 of the level-5
table. The original set is unchanged; this folder exists so the two can be
run and compared without touching it.

## What differs

| | `../hcn_hcl/` | this folder |
|---|---|---|
| chlorine DFT grid | level 5 table, 110 × 770 (pruned) | `CL_ATOM_GRID = (400, 770)` (pruned) |
| other atoms | level 5 table | level 5 table |
| sum-rule repair | on | on, as a check |
| file names | `….pb4d.<date>.py` | `….pb4d.cl400.<date>.py` |
| input directory on the cluster | `dimer-alpha-inputs/cneodft/hcnhcl` | `dimer-alpha-inputs/cneodft/hcnhcl-densegrid` |
| output directory | `dimer-alpha-outputs/cneodft/hcnhcl` | `dimer-alpha-outputs/cneodft/hcnhcl-densegrid` |
| SLURM job names | `hcnhcl-cneo-dft-<tier>` | `hcnhcl-densegrid-cneo-dft-<tier>` |

Everything else, the 45 starting geometries (same file, same PBE
substitutes for the six BLYP and BP86 halogen-bonded starts), the SCF,
optimizer and CPHF settings, the summary printout and the resource tiers,
is identical. `diff` of an input against its `../hcn_hcl/` counterpart
shows only the docstring paragraph, the `CL_ATOM_GRID` setting, its use
in `build_mf`, and two extra grid lines in the printout.

## Why

PySCF's analytic Hessian holds the DFT grid fixed, so its
exchange-correlation part violates the translational sum rule by the
quadrature stiffness of each atom's core density on its own radial grid.
The February 2026 HCN···HCl Hessians (grid level 3, 80 radial shells on
Cl) were off by 0.07 to 0.4 Eh/Bohr² on the chlorine row, which produced
the spurious imaginary bends the original set repairs. The HCN···HBr probe
(`../hcn_hbr/probe/RESULTS.md`) showed that the defect is radial
quadrature error and nothing else: for bromine it fell from 0.73 to
6.7e-4 to 9e-7 Eh/Bohr² for 120, 240 and 400 radial shells, pruning and
angular order had no effect, and the diagonal-block repair reproduced the
converged spectrum to 0.1 cm⁻¹. Chlorine's core is less compact than
bromine's, so 400 shells is generous; it is chosen for uniformity with the
HBr set rather than tuned. Each log prints the per-atom violation before
and after the repair, so the chlorine number shows directly how far the
raw Hessian is from the sum rule. On the original set at level 5 it has
not been measured; the February level-3 values were 0.07 to 0.4, and the
bromine probe gained only a factor of 3 from level 3 to level 5, so
expect 1e-2 to 1e-1 there. Here it should sit near 1e-6.

The variant costs more grid points (roughly 60 % more, the chlorine share
scaling with its radial count) and correspondingly more time in the
grid-dependent parts of the SCF, gradient and Hessian; the tiers have a
margin of more than an order of magnitude for this (see
`../hcn_hbr/README.md`, "Cost").

## Running

```bash
./submit_array.slurm
```

Copy the inputs to `dimer-alpha-inputs/cneodft/hcnhcl-densegrid` first;
the script writes to `dimer-alpha-outputs/cneodft/hcnhcl-densegrid`. The
`cl400` token in every file name keeps the logs distinct from the original
set even if both are ever placed in one directory. Regenerate with
`python make_inputs.py [date]`; `CL_ATOM_GRID = None` in the generator
reproduces the original set apart from the token.
