# (HCl)₂ — conventional DFT re-optimization, Hessian and harmonic frequencies

The all-classical-nuclei counterpart of `../hcl_hcl/`. Thirty self-contained
PySCF drivers: five functionals (PBE, PW91, BP86, BLYP, B97) × three
electronic bases (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ) × two starts for
the hydrogen chloride dimer, conventional Kohn–Sham DFT (`pyscf.dft.RKS`),
no quantum nuclei, 400 radial shells on chlorine as in `../hcl_hcl/`. Each
input

1. runs a tight DFT SCF at the starting geometry,
2. re-optimizes it with geomeTRIC and the analytic DFT gradient,
3. runs a single point at the optimized geometry,
4. computes the analytic DFT Hessian with the hardened CPHF solve,
5. checks the Hessian's translational sum rule, repairs its diagonal
   blocks, and does the harmonic analysis on both the raw and the
   repaired Hessian, and
6. prints a summary: method, functional, electronic basis, nuclear masses,
   the DFT grid actually used (level, radial scheme, partition, pruning,
   radial × angular points per element, total after pruning), the
   optimized geometry (symbols and Cartesians, Angstrom), the key
   distances and angles with a C₂h symmetry measure, the single-point
   energy, the sum-rule violation, the number of imaginary modes raw and
   repaired, all 3N frequencies with the translations/rotations tagged,
   and the vibrational frequencies with those modes dropped.

```bash
python HClHCl.dft-pbe.global-min.aug-cc-pvdz.2026-09-20.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HClHCl.dft-<xc>.<global-min|c2h-saddle>.<electronic basis>.<date>.py
```

No `pb4d` token, since there is no protonic basis.

| Token | Structure | Expected result |
|---|---|---|
| `global-min` | bent C_s Cl–H···Cl–H, one HCl the donor, the acceptor nearly perpendicular; Cl···Cl 3.69 to 3.89 Å, H···Cl 2.39 to 2.60 Å, Cl–H···Cl 172 to 176°, acceptor tilt 92 to 95° in the April DFT runs | minimum, 6 real modes |
| `c2h-saddle` | C₂h structure with both H–Cl bonds tilted 50° from the Cl···Cl axis in a trans arrangement, Cl···Cl 3.85 Å, the donor–acceptor interchange saddle point | saddle point, 1 imaginary mode, if the optimizer keeps the symmetry |

Atom order is H0 (donor H), Cl1 (acceptor Cl), H2 (acceptor H), Cl3 (donor
Cl) in every file, the convention of `../hcl_hcl/` and `../hf_hf/`; H0–Cl3
and Cl1–H2 are the covalent bonds, H0···Cl1 the hydrogen bond.

## Starting geometries

`geometries/hcl_hcl_dft_start_geometries.xyz` holds the 30 starting points
with provenance in each comment line; `make_geometries.py` regenerates it
from the April logs (it reads only the logs whose method token is `dft`).

- The 15 `global-min` starts are the final geometries of the April 2026 DFT
  runs of the same functional and basis. All 15 converged (13 to 16 steps)
  and all are minima (see `april_analysis/`). In the April logs the donor
  was H2 (H2–Cl3 donating to Cl1), as in the CNEO April runs, so the
  extractor swaps the two hydrogens to put the donor at H0 and records the
  swap in the comment line.
- The 15 `c2h-saddle` starts are constructed: Cl···Cl 3.85 Å, tilt 50°,
  r(H–Cl) equal to the acceptor H–Cl length of the April DFT geometry of
  the same functional and basis, exact inversion symmetry to floating-point
  precision, the same construction as in `../hcl_hcl/`. An unconstrained
  minimizer started there has no symmetry-breaking gradient component, so
  it stays on the C₂h surface and converges to the saddle point; the
  printout reports the C₂h asymmetry so the outcome is unambiguous either
  way.

## What the April 2026 DFT runs had, and what changed

The classification in `april_analysis/` (`hcl_hcl_dft_classification.txt`)
found:

1. **All 15 global-minimum runs are minima, raw and repaired.** The
   Hessians finished in 11 to 13 CPHF iterations.
2. **The second start was never defined.** The `lm_start_geom_hydb`
   placeholder holds `DNE` in every field, still with the `F` symbols of the
   (HF)₂ template, so all 15 local-minimum runs died while parsing the
   geometry (two-line logs), as in the CNEO set.
3. **The grid was level 9** (200 radial shells and 1454 angular points on
   every atom), not the level 3 of the CNEO-DFT runs. The chlorine sum-rule
   violation stayed below 1e-4 Eh/Bohr² and the repair moved no vibration
   by more than 0.6 cm⁻¹. The April CNEO Hessians at level 3 (80 shells on
   Cl) carried 0.07 to 0.42 Eh/Bohr² and their raw aVDZ and aVQZ spectra
   showed a spurious 206i to 555i mode.
4. **The optimizer settings** were the custom set of all April runs
   (grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3 Å, ΔE 1e-9, 1500-step
   cap) with default SCF tolerances and no grid response. The starting
   coordinates were genuinely in Bohr.

## CNEO-DFT minus DFT (April runs, repaired Hessians, H 1.008 u on both sides)

The April CNEO-DFT (HCl)₂ runs had only the acceptor proton quantum, the
reverse of (HF)₂. Against the DFT runs, CNEO-DFT lowers the acceptor H–Cl
stretch by 103 to 119 cm⁻¹ for PBE, PW91 and BP86 and the donor stretch, a
classical proton in both, by 8 to 10 cm⁻¹. For BLYP and B97 the lowered
acceptor stretch crosses the donor stretch, so rank matching splits the
shift over the two stretches (16 and 93 cm⁻¹ for BLYP/aVDZ, 53 and 70 for
B97/aVQZ); their sum, 108 to 125 cm⁻¹, is the quantity to compare. The
four intermolecular modes change by only −2 to +8 cm⁻¹, H···Cl shortens by
0.010 to 0.015 Å and the acceptor H–Cl lengthens by 0.023 to 0.027 Å: a
quantum acceptor proton hardly touches the hydrogen bond, as in HF/HCl.
The new CNEO inputs in `../hcl_hcl/` have both protons quantum. Part 4 of
the classification lists every pair.

## Settings

Those of `../hcl_hcl/`, applied to `dft.RKS`: SCF `conv_tol` 1e-11 /
`conv_tol_grad` 1e-6 / 200 cycles, grid level 5 (H 70 × 590) with 400
radial shells × 770 angular points on chlorine
(`ATOM_GRIDS = {'Cl': (400, 770)}`, the level-5 table would give 110), grid
response in the gradient, geomeTRIC GAU_TIGHT with `convergence_energy`
1e-8 and `subfrctor=2`, 200-step cap, `STOP_IF_OPT_FAILS`, `conv_tol_cpscf`
1e-9, 300 Krylov iterations with the level-shift retry, the sum-rule repair
with the raw spectrum printed alongside, and a dissociation warning at
Cl···Cl above 6.5 Å. `ATOM_GRIDS = {}` would fall back to the plain level-5
table.

**Masses.** PySCF's isotope-averaged atomic masses (H 1.008, Cl 35.45 u),
the standard convention for conventional DFT and the one the April logs
used; the summary prints them. The CNEO mass convention of the CNEO-DFT
inputs (nuclear mass 1.007276 u for the quantum protons) is not used in
these files; the 0.07 % difference in the hydrogen mass changes a frequency
by at most 0.04 %, about 1 cm⁻¹ at 3000 cm⁻¹.

**Why the dense chlorine grid and the repair are kept for conventional
DFT.** PySCF's RKS Hessian holds the grid fixed for GGA and hybrid-GGA
functionals exactly as the NEO Hessian does, so the chlorine diagonal-block
defect appears here too; at level 3 it produced the spurious imaginary
modes of the CNEO runs. The HCN···HBr probe showed that 400 radial shells
remove it and that the repair reproduces the converged spectrum to
0.1 cm⁻¹. Every log prints the per-atom violation before and after the
repair; both chlorine numbers should sit near 1e-6.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (6 for
both structures); PySCF's projected analysis as a one-line cross-check.
Both tables use the repaired Hessian when the repair is on.

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/dft/hclhcl`, writing all logs
and results to `.../dimer-alpha-outputs/dft/hclhcl`; change the `dft`
directory token in the script if the DFT runs live elsewhere. The tiers are
those of the CNEO-DFT sets and are generous by more than an order of
magnitude for a four-atom conventional-DFT calculation. To run the minima
alone, leave the `c2h-saddle` files out of the input directory. Regenerate
the inputs with `python make_inputs.py [date]` after editing the template
in `make_inputs.py`.

## april_analysis/

`analyze_hxhx_dft.py` (with its parser `reanalyze_hessians2.py`), the same
generic (HX)₂ classifier as in `../hf_hf_dft/`; `hcl_hcl_dft_classification.txt`
is its output for the 30 non-5Z DFT logs (15 with a Hessian), with
CNEO-DFT minus DFT shifts in Part 4 and the quantum proton of each CNEO run
noted. Run it as `python analyze_hxhx_dft.py <DFT log dir> [<CNEO-DFT log dir>]`.
