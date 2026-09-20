# (HF)₂ — conventional DFT re-optimization, Hessian and harmonic frequencies

The all-classical-nuclei counterpart of `../hf_hf/`. Thirty self-contained
PySCF drivers: five functionals (PBE, PW91, BP86, BLYP, B97) × three
electronic bases (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ) × two starts for
the hydrogen fluoride dimer, conventional Kohn–Sham DFT (`pyscf.dft.RKS`),
no quantum nuclei, 400 radial shells on fluorine as in `../hf_hf/`. Each
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
python HFHF.dft-pbe.global-min.aug-cc-pvdz.2026-09-20.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HFHF.dft-<xc>.<global-min|c2h-saddle>.<electronic basis>.<date>.py
```

No `pb4d` token, since there is no protonic basis.

| Token | Structure | Expected result |
|---|---|---|
| `global-min` | bent C_s F–H···F–H, one HF the donor, the acceptor tilted; F···F 2.71 to 2.77 Å, H···F 1.77 to 1.86 Å, F–H···F 169 to 171°, acceptor tilt 109 to 114° in the April DFT runs | minimum, 6 real modes |
| `c2h-saddle` | C₂h structure with both H–F bonds tilted 55° from the F···F axis in a trans arrangement, F···F 2.78 Å, the donor–acceptor interchange saddle point | saddle point, 1 imaginary mode, if the optimizer keeps the symmetry |

Atom order is H0 (donor H), F1 (acceptor F), H2 (acceptor H), F3 (donor F)
in every file, the convention of `../hf_hf/`; H0–F3 and F1–H2 are the
covalent bonds, H0···F1 the hydrogen bond.

## Starting geometries

`geometries/hf_hf_dft_start_geometries.xyz` holds the 30 starting points
with provenance in each comment line; `make_geometries.py` regenerates it
from the April logs (it reads only the logs whose method token is `dft`).

- The 15 `global-min` starts are the final geometries of the April 2026 DFT
  runs of the same functional and basis. All 15 converged, in 30 to 42
  steps: the April inputs gave the starting coordinates in Ångström but
  declared them as Bohr, so every run began from a structure compressed by
  a factor of 1.89 and spent most of its steps expanding. All 15 are minima
  (see `april_analysis/`). The donor was H0 in every April geometry, so no
  reordering was needed; the extractor checks this from the distances.
- The 15 `c2h-saddle` starts are constructed: F···F 2.78 Å, tilt 55°,
  r(H–F) equal to the acceptor H–F length of the April DFT geometry of the
  same functional and basis, exact inversion symmetry to floating-point
  precision, the same construction as in `../hf_hf/`. An unconstrained
  minimizer started there has no symmetry-breaking gradient component, so
  it stays on the C₂h surface and converges to the saddle point; the
  printout reports the C₂h asymmetry so the outcome is unambiguous either
  way.

## What the April 2026 DFT runs had, and what changed

The classification in `april_analysis/` (`hf_hf_dft_classification.txt`)
found:

1. **All 15 global-minimum runs are minima, raw and repaired.** The
   Hessians finished in 11 CPHF iterations each.
2. **The second start was never defined.** The `lm_start_geom_hydb`
   placeholder holds `DNE` in every field, so all 15 local-minimum runs
   died while parsing the geometry (two-line logs), as in the CNEO set.
3. **The grid was level 9** (200 radial shells and 1454 angular points on
   every atom), not the level 3 of the CNEO-DFT runs. The fluorine sum-rule
   violation was 1e-6 to 2e-4 Eh/Bohr² and the hydrogen violation up to
   2.4e-4; the repair moved no vibration by more than 6 cm⁻¹ (median 0.7).
4. **The optimizer settings** were the custom set of all April runs
   (grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3 Å, ΔE 1e-9, 1500-step
   cap) with default SCF tolerances and no grid response.

## CNEO-DFT minus DFT (April runs, repaired Hessians, H 1.008 u on both sides)

The April CNEO-DFT (HF)₂ runs had only the donor proton quantum. Against
the DFT runs, CNEO-DFT lowers the donor H–F stretch by 206 to 240 cm⁻¹ and
the acceptor H–F stretch, whose proton was classical in both, by 6 to
7 cm⁻¹; the four intermolecular modes rise by 6 to 38 cm⁻¹, the donor H–F
bond lengthens by 0.023 to 0.025 Å, H···F shortens by 0.034 to 0.046 Å
and F···F by 0.010 to 0.022 Å: the quantum donor proton strengthens the
hydrogen bond. Part 4 of the classification lists every pair. The new
CNEO inputs in `../hf_hf/` have both protons quantum, so the acceptor
stretch will shift as well when they run.

## Settings

Those of `../hf_hf/`, applied to `dft.RKS`: SCF `conv_tol` 1e-11 /
`conv_tol_grad` 1e-6 / 200 cycles, grid level 5 (H 70 × 590) with 400
radial shells × 770 angular points on fluorine
(`ATOM_GRIDS = {'F': (400, 770)}`, the level-5 table would give 105), grid
response in the gradient, geomeTRIC GAU_TIGHT with `convergence_energy`
1e-8 and `subfrctor=2`, 200-step cap, `STOP_IF_OPT_FAILS`, `conv_tol_cpscf`
1e-9, 300 Krylov iterations with the level-shift retry, the sum-rule repair
with the raw spectrum printed alongside. `ATOM_GRIDS = {}` would fall back
to the plain level-5 table.

**Masses.** PySCF's isotope-averaged atomic masses (H 1.008, F 18.998 u),
the standard convention for conventional DFT and the one the April logs
used; the summary prints them. The CNEO mass convention of the CNEO-DFT
inputs (nuclear mass 1.007276 u for the quantum protons) is not used in
these files; the 0.07 % difference in the hydrogen mass changes a frequency
by at most 0.04 %, about 1 cm⁻¹ at 3000 cm⁻¹.

**Why the dense fluorine grid and the repair are kept for conventional
DFT.** PySCF's RKS Hessian holds the grid fixed for GGA and hybrid-GGA
functionals exactly as the NEO Hessian does, so the fluorine diagonal-block
defect appears here too; at level 3 it moved the two lowest CNEO modes by
12 to 22 cm⁻¹ at aug-cc-pVQZ. The probe showed that 400 radial shells
remove it and that the repair reproduces the converged spectrum to
0.1 cm⁻¹. Every log prints the per-atom violation before and after the
repair.

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
`.../dimers/pyscf-master/dimer-alpha-inputs/dft/hfhf`, writing all logs
and results to `.../dimer-alpha-outputs/dft/hfhf`; change the `dft`
directory token in the script if the DFT runs live elsewhere. The tiers are
those of the CNEO-DFT sets and are generous by more than an order of
magnitude for a four-atom conventional-DFT calculation. To run the minima
alone, leave the `c2h-saddle` files out of the input directory. Regenerate
the inputs with `python make_inputs.py [date]` after editing the template
in `make_inputs.py`.

## april_analysis/

`analyze_hxhx_dft.py` (with its parser `reanalyze_hessians2.py`) serves
every (HX)₂ set: it reads the halogen from the logs, re-diagonalizes the
Hessians printed in the April DFT logs with the six-mode projection, reads
the covalent partners and the donor from the distances, checks the sum rule
per atom, applies the repair and classifies each run on the raw and on the
repaired Hessian; given the CNEO-DFT log directory as a second argument it
also prints CNEO-DFT minus DFT frequency and geometry shifts with the
quantum proton of each CNEO run noted. `hf_hf_dft_classification.txt` is
its output for the 30 non-5Z DFT logs (15 with a Hessian). Run it as
`python analyze_hxhx_dft.py <DFT log dir> [<CNEO-DFT log dir>]`.
