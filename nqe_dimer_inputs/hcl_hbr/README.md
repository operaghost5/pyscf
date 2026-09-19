# HCl/HBr — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Thirty self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × two isomers of the hydrogen chloride / hydrogen bromide
dimer, **both** hydrogen nuclei quantum with the PB4-D protonic basis, no
electron-proton correlation functional, both halogens all-electron, 400
radial shells on both chlorine and bromine. Each input

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
   Cartesians, Angstrom), the key distances and angles with the detected
   donor, the single-point energy, the sum-rule violation, the number of
   imaginary modes raw and repaired, all 3N frequencies with the
   translations/rotations tagged, and the vibrational frequencies with
   those modes dropped.

```bash
python HClHBr.cneo-dft-pbe.hcl-donor.aug-cc-pvdz.pb4d.2026-09-19.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HClHBr.cneo-dft-<xc>.<hcl-donor|hbr-donor>.<electronic basis>.pb4d.<date>.py
```

| Token | Structure | April 2026 result (repaired Hessian) |
|---|---|---|
| `hcl-donor` | Cl–H···Br–H, HCl the donor, HBr roughly perpendicular to the hydrogen bond; H···Br 2.46 to 2.66 Å, tilt 90 to 93° | minimum for every functional and basis that finished |
| `hbr-donor` | Br–H···Cl–H, HBr the donor, HCl roughly perpendicular; H···Cl 2.38 to 2.61 Å, tilt 91 to 95° | minimum for every functional and basis |

The April inputs called these `global_min` and `local_min`; the donor-based
names are used because the April energies cannot rank the isomers (see
below). Atom order is H0, Cl1, H2, Br3 in every file (H0–Cl1 is HCl, H2–Br3
is HBr); atoms 0 and 2 are the quantum protons.

## Starting geometries

`geometries/hcl_hbr_start_geometries.xyz` holds the 30 starting points, the
final geometries of the April 2026 runs of the same functional, basis and
isomer, with the run's energy, convergence status and quantum nucleus in
each comment line. 29 of the April optimizations converged (10 to 18
steps, B97 37 to 57). The thirtieth, B97/aug-cc-pVQZ `hcl-donor`, was
killed after 438 steps while oscillating at a gradient norm of 4e-5
Eh/Bohr, the net-force artefact that `subfrctor=2` removes; its start is
the geometry of that last completed step, which the comment line records.
`make_geometries.py` regenerates the file from the April logs, decides
from the distances which hydrogen belongs to which molecule, and falls
back to the last optimizer step when a log has no final geometry.

## What the April 2026 runs had, and what changed

The classification in `april_analysis/` found:

1. **Only the donor proton was quantum, and a different one in each
   isomer.** `quantum_nuc=[0]` (the HCl proton) for the HCl-donor runs,
   `quantum_nuc=[2]` (the HBr proton) for the HBr-donor runs. The two sets
   of energies carry different nuclear quantum contributions, so their
   differences (±0.25 kcal/mol) say nothing about which isomer is lower.
   These inputs use `QUANTUM_NUC = ['H']`, both protons, so the
   re-optimized energies rank the isomers directly.
2. **The bromine sum-rule defect dominated the raw spectra.** At grid level
   3 the bromine row was off by 1.4 to 2.7 Eh/Bohr² (aVDZ), 9.5 to 18
   (aVTZ) and 2.7 to 3.3 (aVQZ), chlorine by 0.02 to 0.42. A displaced
   translation leaked into the printed vibrations: a spurious 1000 to
   1390 cm⁻¹ mode at TZ, a spurious ~700i mode at QZ that made all nine
   QZ runs look like first-order saddle points, and the lowest real mode
   pushed out of the list at DZ. After the diagonal-block repair all 29
   spectra have six real modes and the lowest mode agrees across bases to
   about 1 cm⁻¹ per functional. Both isomers are minima everywhere.
3. **Grid level 3 and the loose optimizer settings**, as in the other
   April sets. The starting coordinates were genuinely in Bohr.

## Settings

Those of the other sets (see `../hcn_hf/README.md` and
`../hcn_hbr/README.md`): GAU_TIGHT criteria with `subfrctor=2`, SCF
`conv_tol` 1e-11 / `conv_tol_grad` 1e-6, grid level 5 with grid response
on, 200-step cap, `STOP_IF_OPT_FAILS`, CPHF with 300 Krylov iterations and
a level-shift retry, the sum-rule repair with the raw spectrum printed
alongside.

**Dense halogen radial grids.** `ATOM_GRIDS = {'Cl': (400, 770), 'Br':
(400, 770)}` puts 400 radial shells on both halogens, where the level-5
table gives Cl 110 and Br 120; hydrogen keeps the table. The HCN···HBr
probe (`../hcn_hbr/probe/RESULTS.md`) showed that the halogen defect is
radial quadrature error that 400 shells remove and that the repair
reproduces the converged spectrum to 0.1 cm⁻¹. Every log prints the
per-atom violation before and after the repair; both halogen numbers
should sit near 1e-6. `ATOM_GRIDS = {}` restores the plain level-5 grid.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (6 for
both isomers); PySCF's projected analysis as a one-line cross-check. Both
tables use the repaired Hessian when the repair is on. Masses are the CNEO
nuclear masses (`mol.mass`).

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/cneodft/hclhbr`, writing all
logs and results to `.../dimer-alpha-outputs/cneodft/hclhbr`. The tiers are
those of the other sets; with four atoms and at most 260 electronic AOs
they are generous by more than an order of magnitude. Regenerate the
inputs with `python make_inputs.py [date]` after editing the template in
`make_inputs.py`.

## april_analysis/

`analyze_hclbr.py` (with its parser `reanalyze_hessians2.py`)
re-diagonalizes the Hessians printed in the April logs with the six-mode
projection, reads the bonding from the distances, checks the sum rule per
atom, applies the repair and classifies each run on the raw and on the
repaired Hessian; `hcl_hbr_classification.txt` is its output for the 31
non-5Z logs (29 with a Hessian). Run it as
`python analyze_hclbr.py <April log dir>`.
