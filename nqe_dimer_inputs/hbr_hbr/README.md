# (HBr)₂ — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Thirty self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × two starts for the hydrogen bromide dimer, **both**
hydrogen nuclei quantum with the PB4-D protonic basis, no electron-proton
correlation functional, bromine all-electron, 400 radial shells on
bromine. Each input

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
python HBrHBr.cneo-dft-pbe.global-min.aug-cc-pvdz.pb4d.2026-09-19.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HBrHBr.cneo-dft-<xc>.<global-min|c2h-saddle>.<electronic basis>.pb4d.<date>.py
```

| Token | Structure | Expected result |
|---|---|---|
| `global-min` | bent C_s Br–H···Br–H, one HBr the donor, the acceptor nearly perpendicular; Br···Br 3.99 to 4.23 Å, H···Br 2.54 to 2.80 Å, Br–H···Br 174 to 177°, acceptor tilt 89 to 92° in the April runs | minimum, 6 real modes |
| `c2h-saddle` | C₂h structure with both H–Br bonds tilted 50° from the Br···Br axis in a trans arrangement, Br···Br 4.15 Å, the donor–acceptor interchange saddle point | saddle point, 1 imaginary mode, if the optimizer keeps the symmetry |

Atom order is H0 (donor H), Br1 (acceptor Br), H2 (acceptor H), Br3 (donor
Br) in every file, the same convention as `../hf_hf/` and `../hcl_hcl/`;
H0–Br3 and Br1–H2 are the covalent bonds, H0···Br1 the hydrogen bond. Atoms 0 and 2 are the
quantum protons.

## Starting geometries

`geometries/hbr_hbr_start_geometries.xyz` holds the 30 starting points with
provenance in each comment line; `make_geometries.py` regenerates it from
the April logs.

- The 15 `global-min` starts are the final geometries of the April 2026
  runs of the same functional and basis. All 15 converged (11 to 18
  steps; B97 26, 90 and 50) and all are minima on the repaired Hessian
  (see `april_analysis/`). In the April logs the donor was H2 (H2–Br3
  donating to Br1), so the two hydrogens are swapped to put the donor at
  H0; the extractor decides the donor from the distances and records the
  swap.
- The 15 `c2h-saddle` starts are constructed: Br···Br 4.15 Å, tilt 50°,
  r(H–Br) equal to the acceptor H–Br length of the April geometry of the
  same functional and basis, exact inversion symmetry to floating-point
  precision. An unconstrained minimizer started there has no
  symmetry-breaking gradient component, so it stays on the C₂h surface and
  converges to the saddle point, whose Hessian then has one imaginary
  mode. If numerical noise breaks the symmetry the run slides down to the
  bent minimum; the printout reports the C₂h asymmetry (largest difference
  between the two H–Br or the two H···Br distances), so the outcome is
  unambiguous either way.

## What the April 2026 runs had, and what changed

The classification in `april_analysis/` found:

1. **Only the acceptor proton was quantum.** `quantum_nuc=[0]` selects H0,
   but the distances show H2 donating the hydrogen bond in every run, as in
   (HCl)₂ and unlike (HF)₂, where the donor was quantum. These inputs use
   `QUANTUM_NUC = ['H']`, both protons.
2. **The second start was never defined.** The `lm_start_geom_hydb`
   placeholder holds `DNE` in every field, still with the `F` symbols of the
   (HF)₂ template, so all 20 local-minimum runs died while parsing the
   geometry (two-line logs).
3. **The bromine sum-rule defect made every raw spectrum wrong.** At grid
   level 3 the bromine rows were off by 1.4 to 2.7 Eh/Bohr² (aVDZ), 9.5 to
   18 (aVTZ) and 2.7 to 3.3 (aVQZ). A displaced translation appeared as a
   spurious 690 to 945 cm⁻¹ mode at DZ (replacing the lowest real mode), as
   a spurious 1766 to 2327 cm⁻¹ mode among the H–Br stretches at TZ (with
   the three lowest modes shifted by 100 to 200 cm⁻¹), and as a spurious
   934i to 1045i mode at QZ that made all five QZ runs look like saddle
   points (B97 second order, with a 29i companion). After the diagonal-block
   repair all 15 spectra have six real modes and the lowest mode agrees
   across bases to about 2 cm⁻¹ per functional (51.9, 50.6, 50.1 cm⁻¹ for
   PBE).
4. **Grid level 3 and the loose optimizer settings**, as in the other
   April sets. The starting coordinates were genuinely in Bohr.

## Settings

Those of the other sets (see `../hcn_hf/README.md` and
`../hcn_hbr/README.md`): GAU_TIGHT criteria with `subfrctor=2`, SCF
`conv_tol` 1e-11 / `conv_tol_grad` 1e-6, grid level 5 with grid response
on, 200-step cap, `STOP_IF_OPT_FAILS`, CPHF with 300 Krylov iterations and
a level-shift retry, the sum-rule repair with the raw spectrum printed
alongside.

**Dense bromine radial grid.** `BR_ATOM_GRID = (400, 770)` puts 400 radial
shells on bromine where the level-5 table gives 120; hydrogen keeps the
table. The HCN···HBr probe (`../hcn_hbr/probe/RESULTS.md`) measured exactly
this: 0.73, 6.7e-4 and 9e-7 Eh/Bohr² for 120, 240 and 400 shells, and the
repair reproducing the converged spectrum to 0.1 cm⁻¹. Every log prints
the per-atom violation before and after the repair; both bromine numbers
should sit near 1e-6. `BR_ATOM_GRID = None` restores the plain level-5
grid.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (6 for
both structures); PySCF's projected analysis as a one-line cross-check.
Both tables use the repaired Hessian when the repair is on. Masses are the
CNEO nuclear masses (`mol.mass`).

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/cneodft/hbrhbr`, writing all
logs and results to `.../dimer-alpha-outputs/cneodft/hbrhbr`. The tiers are
those of the other sets; with four atoms and at most 260 electronic AOs
they are generous by more than an order of magnitude. To run the minima
alone, leave the `c2h-saddle` files out of the input directory. Regenerate
the inputs with `python make_inputs.py [date]` after editing the template
in `make_inputs.py`.

## april_analysis/

`analyze_hbrbr.py` (with its parser `reanalyze_hessians2.py`)
re-diagonalizes the Hessians printed in the April logs with the six-mode
projection, reads the donor from the distances, checks the bromine sum
rule, applies the repair and classifies each run on the raw and on the
repaired Hessian; `hbr_hbr_classification.txt` is its output for the 30
non-5Z logs. Run it as `python analyze_hbrbr.py <April log dir>`.
