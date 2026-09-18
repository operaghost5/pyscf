# HF/HBr — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Thirty self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × two isomers of the hydrogen fluoride / hydrogen bromide
dimer, **both** hydrogen nuclei quantum with the PB4-D protonic basis, no
electron-proton correlation functional, bromine all-electron, 400 radial
shells on both bromine and fluorine. Each input

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
python HFHBr.cneo-dft-pbe.hf-donor.aug-cc-pvdz.pb4d.2026-09-18.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HFHBr.cneo-dft-<xc>.<hf-donor|hbr-donor>.<electronic basis>.pb4d.<date>.py
```

| Token | Structure | April 2026 result (repaired Hessian) |
|---|---|---|
| `hf-donor` | F–H···Br–H, HF the donor, HBr roughly perpendicular to the hydrogen bond; H···Br 2.32 to 2.46 Å, tilt 91 to 92° | minimum for every functional and basis |
| `hbr-donor` | Br–H···F–H, HBr the donor, HF tilted; H···F 1.96 to 2.14 Å, tilt 108 to 117° | minimum for every functional and basis |

The April inputs called these `global_min` and `local_min`. The labels are
not used here because the April energies cannot rank the isomers (see
below); the two are named by the donor instead. Atom order is H0, F1, H2,
Br3 in every file (H0–F1 is HF, H2–Br3 is HBr); atoms 0 and 2 are the
quantum protons.

## Starting geometries

`geometries/hf_hbr_start_geometries.xyz` holds the 30 starting points, the
final geometries of the April 2026 runs of the same functional, basis and
isomer, with the run's energy, convergence status and quantum nucleus in
each comment line. All 30 April optimizations converged (12 to 28 steps,
B97 up to 191). `make_geometries.py` regenerates the file from the April
logs.

## What the April 2026 runs had, and what changed

The classification in `april_analysis/` found three things:

1. **Only the donor proton was quantum, and a different one in each
   isomer.** `quantum_nuc=[0]` (the HF proton) for the HF-donor runs,
   `quantum_nuc=[2]` (the HBr proton) for the HBr-donor runs. The two sets
   of energies therefore carry different nuclear quantum contributions of
   several kcal/mol each, and their differences (±0.1 kcal/mol, −0.5 for
   B97) say nothing about which isomer is lower. These inputs use
   `QUANTUM_NUC = ['H']`, both protons, so the re-optimized energies rank
   the isomers directly.
2. **The bromine sum-rule defect dominated the raw spectra.** At grid level
   3 the bromine row was off by 1.4 to 2.7 Eh/Bohr² (aVDZ), 9.5 to 18
   (aVTZ) and 2.7 to 3.3 (aVQZ). A displaced translation leaked into the
   printed vibrations: the lowest real mode near 110 cm⁻¹ was replaced by
   a spurious one near 410 at DZ, a spurious 800 to 1100 cm⁻¹ mode
   appeared at TZ, and a spurious ~400i mode at QZ made all ten QZ runs
   look like first-order saddle points. After the diagonal-block repair
   all 30 spectra have six real modes and agree across bases to 10 to
   20 cm⁻¹ per functional. Both isomers are minima everywhere.
3. **Grid level 3 and the loose optimizer settings**, as in the other
   April sets. The starting coordinates were genuinely in Bohr here.

## Settings

Those of the other sets (see `../hcn_hf/README.md` and
`../hcn_hbr/README.md`): GAU_TIGHT criteria with `subfrctor=2`, SCF
`conv_tol` 1e-11 / `conv_tol_grad` 1e-6, grid level 5 with grid response
on, 200-step cap, `STOP_IF_OPT_FAILS`, CPHF with 300 Krylov iterations and
a level-shift retry, the sum-rule repair with the raw spectrum printed
alongside.

**Dense halogen radial grids.** `ATOM_GRIDS = {'F': (400, 770), 'Br':
(400, 770)}` puts 400 radial shells on both halogens, where the level-5
table gives F 105 and Br 120; hydrogen keeps the table. The HCN···HBr probe
(`../hcn_hbr/probe/RESULTS.md`) showed that the halogen defect is radial
quadrature error that 400 shells remove and that the repair reproduces
the converged spectrum to 0.1 cm⁻¹. Every log prints the per-atom
violation before and after the repair; both halogen numbers should sit
near 1e-6. `ATOM_GRIDS = {}` restores the plain level-5 grid.

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
`.../dimers/pyscf-master/dimer-alpha-inputs/cneodft/hfhbr`, writing all
logs and results to `.../dimer-alpha-outputs/cneodft/hfhbr`. The tiers are
those of the other sets; with four atoms and at most 256 electronic AOs
they are generous by more than an order of magnitude. Regenerate the
inputs with `python make_inputs.py [date]` after editing the template in
`make_inputs.py`.

## april_analysis/

`analyze_hfbr.py` (with its parser `reanalyze_hessians2.py`) re-diagonalizes
the Hessians printed in the April logs with the six-mode projection, checks
the sum rule per atom, applies the repair and classifies each run on the
raw and on the repaired Hessian; `hf_hbr_classification.txt` is its output
for the 30 non-5Z logs. Run it as `python analyze_hfbr.py <April log dir>`.
