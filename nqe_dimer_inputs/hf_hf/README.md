# (HF)₂ — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Thirty self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × two starts for the hydrogen fluoride dimer, **both**
hydrogen nuclei quantum with the PB4-D protonic basis, no electron-proton
correlation functional, 400 radial shells on fluorine. Each input

1. runs a tight CNEO-DFT SCF at the starting geometry,
2. re-optimizes it with geomeTRIC and the analytic CNEO-DFT gradient,
3. runs a single point at the optimized geometry,
4. computes the analytic CNEO-DFT Hessian with the hardened CPHF solve,
5. checks the Hessian's translational sum rule, repairs its diagonal
   blocks, and does the harmonic analysis on both the raw and the
   repaired Hessian, and
6. prints a summary: method, functional, quantum nuclei, EPC, electronic
   and nuclear basis, nuclear masses, the DFT grid actually used (level,
   radial scheme, partition, pruning, radial × angular points per element
   and the total after pruning), the optimized geometry (symbols and
   Cartesians, Angstrom), the key distances and angles with a C₂h symmetry
   measure, the single-point energy, the sum-rule violation, the number of
   imaginary modes raw and repaired, all 3N frequencies with the
   translations/rotations tagged, and the vibrational frequencies with
   those modes dropped.

```bash
python HFHF.cneo-dft-pbe.global-min.aug-cc-pvdz.pb4d.2026-09-18.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HFHF.cneo-dft-<xc>.<global-min|c2h-saddle>.<electronic basis>.pb4d.<date>.py
```

| Token | Structure | Expected result |
|---|---|---|
| `global-min` | bent C_s F–H···F–H, one HF the donor; F···F 2.69 to 2.75 Å, H···F 1.73 to 1.81 Å, F–H···F 169 to 170°, acceptor tilt 109 to 114° in the April runs | minimum, 6 real modes |
| `c2h-saddle` | C₂h structure with both H–F bonds tilted 55° from the F···F axis in a trans arrangement, the donor–acceptor interchange saddle point (about 1 kcal/mol above the minimum on the electronic surface) | saddle point, 1 imaginary mode, if the optimizer keeps the symmetry |

Atom order is H0 (donor H), F1 (acceptor F), H2 (acceptor H), F3 (donor
F) in every file; H0–F3 and F1–H2 are the covalent bonds, H0···F1 the
hydrogen bond. Atoms 0 and 2 are the quantum protons.

## Starting geometries

`geometries/hf_hf_start_geometries.xyz` holds the 30 starting points with
provenance in each comment line; `make_geometries.py` regenerates it from
the April logs.

- The 15 `global-min` starts are the final geometries of the April 2026
  runs of the same functional and basis. All 15 converged (36 to 103
  steps) and all are minima (see `april_analysis/`).
- The 15 `c2h-saddle` starts are constructed: F···F 2.78 Å, tilt 55°,
  r(H–F) equal to the acceptor H–F length of the April geometry of the same
  functional and basis. They have exact inversion symmetry to
  floating-point precision. An unconstrained minimizer started there has
  no symmetry-breaking gradient component, so it stays on the C₂h surface
  and converges to the saddle point, whose Hessian then has one imaginary
  mode. If numerical noise breaks the symmetry along the way the run
  slides down to the bent minimum instead; the printout reports the C₂h
  asymmetry (largest difference between the two H–F or the two H···F
  distances), so the outcome is unambiguous either way. There is no
  second minimum on the (HF)₂ surface to look for.

## What the April 2026 runs had, and what changed

The April inputs (`hydrogen_fluoride_hydrogen_fluoride_dimer.*`) differ
from these in four ways that the classification in `april_analysis/`
uncovered:

1. **Only the donor proton was quantum.** `quantum_nuc=[0]` selects atom 0
   alone; every SCF cycle in the logs carries a single nuclear component.
   The acceptor H–F stretch (3850 to 4070 cm⁻¹) is therefore a classical
   harmonic value while the donor stretch (3490 to 3740 cm⁻¹) carries the
   CNEO lowering. These inputs use `QUANTUM_NUC = ['H']`, both protons, as
   in the HCN···HX sets.
2. **The start was an Ångström geometry read as Bohr.** `unit='Bohr'` with
   H–F 0.92 and F···F 2.73 compressed every start by a factor of 1.89
   (initial SCF energies 1.9 Eh above the final ones). The optimizer
   recovered in all 15 cases. These inputs are in Ångström.
3. **The second start was never defined.** The `lm_start_geom_hydb`
   placeholder held `DNE` in every field, so all 20 local-minimum runs died
   while parsing the geometry, before the first SCF (two-line logs).
4. **Grid level 3.** The fluorine sum-rule violation was 1e-3 Eh/Bohr² at
   aug-cc-pVDZ/TZ and 4e-3 to 9e-3 at aug-cc-pVQZ; at QZ the
   diagonal-block repair moved the two lowest intermolecular modes by 12
   to 22 cm⁻¹ and brought them into agreement with the smaller bases
   (lowest mode within about 3 cm⁻¹ across DZ/TZ/QZ for each functional).
   The classification was unaffected: all 15 are minima with or without
   the repair.

## Settings

Those of the HCN···HF, HCl and HBr sets (see `../hcn_hf/README.md` and
`../hcn_hbr/README.md`): GAU_TIGHT criteria with `subfrctor=2`, SCF
`conv_tol` 1e-11 / `conv_tol_grad` 1e-6, grid level 5 with grid response
on, 200-step cap, `STOP_IF_OPT_FAILS`, CPHF with 300 Krylov iterations and
a level-shift retry, the sum-rule repair with the raw spectrum printed
alongside.

**Dense fluorine radial grid.** `F_ATOM_GRID = (400, 770)` puts 400 radial
shells on fluorine where the level-5 table gives 105; hydrogen keeps the
table. The HCN···HBr probe (`../hcn_hbr/probe/RESULTS.md`) showed that the
halogen sum-rule defect of the fixed-grid Hessian is radial quadrature
error that 400 shells remove, and that the repair reproduces the converged
spectrum to 0.1 cm⁻¹. For a four-atom system the extra grid points cost
little. Every log prints the per-atom violation before and after the
repair; the fluorine numbers should sit near 1e-6. `F_ATOM_GRID = None`
restores the plain level-5 grid.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (6 for
both structures here); PySCF's projected analysis as a one-line
cross-check. Both tables use the repaired Hessian when the repair is on.
Masses are the CNEO nuclear masses (`mol.mass`).

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/cneodft/hfhf`, writing all logs
and results to `.../dimer-alpha-outputs/cneodft/hfhf`. The tiers are those
of the other sets; with four atoms and at most 252 electronic AOs they
are generous by more than an order of magnitude. To run the minima alone,
leave the `c2h-saddle` files out of the input directory. Regenerate the
inputs with `python make_inputs.py [date]` after editing the template in
`make_inputs.py`.

## april_analysis/

`analyze_hfhf.py` (with its parser `reanalyze_hessians2.py`) re-diagonalizes
the Hessians printed in the April logs with the correct six-mode projection,
checks the fluorine sum rule, applies the repair and classifies each run;
`hf_hf_classification.txt` is its output for the 30 non-5Z logs. Run it as
`python analyze_hfhf.py <April log dir>`.
