# HF/HCl — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Thirty self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × two isomers of the hydrogen fluoride / hydrogen chloride
dimer, **both** hydrogen nuclei quantum with the PB4-D protonic basis, no
electron-proton correlation functional, 400 radial shells on both chlorine
and fluorine. Each input

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
python HFHCl.cneo-dft-pbe.hf-donor.aug-cc-pvdz.pb4d.2026-09-19.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HFHCl.cneo-dft-<xc>.<hf-donor|hcl-donor>.<electronic basis>.pb4d.<date>.py
```

| Token | Structure | April 2026 result (repaired Hessian) |
|---|---|---|
| `hf-donor` | F–H···Cl–H, HF the donor, HCl roughly perpendicular to the hydrogen bond; H···Cl 2.21 to 2.34 Å, tilt 92 to 94° | minimum for every functional and basis; lower than `hcl-donor` by 0.07 to 0.38 kcal/mol with the HCl proton quantum and the HF proton classical |
| `hcl-donor` | Cl–H···F–H, HCl the donor, HF tilted; H···F 1.90 to 2.05 Å, tilt 109 to 116° | minimum for every functional and basis |

The April inputs called these `global_min` and `local_min`; the donor-based
names are used here as for HF/HBr. Atom order is H0, F1, H2, Cl3 in every
file (H0–F1 is HF, H2–Cl3 is HCl); atoms 0 and 2 are the quantum protons.

## Starting geometries

`geometries/hf_hcl_start_geometries.xyz` holds the 30 starting points, the
final geometries of the April 2026 runs of the same functional, basis and
isomer, with the run's energy, convergence status and quantum nucleus in
each comment line. All 30 April optimizations converged (11 to 33 steps,
B97 up to 71, and 295 for B97/aug-cc-pVDZ `hf-donor`).

The April `global_min` (HF-donor) inputs had the hydrogens the other way
round, H0 on chlorine and H2 on fluorine; `make_geometries.py` decides from
the distances which hydrogen belongs to which molecule and swaps the two
hydrogens where needed, so that both isomers share one atom order. The
comment line records the swap. Every input also checks at the end that
H0–F1 and H2–Cl3 are still the covalent bonds.

## What the April 2026 runs had, and what changed

The classification in `april_analysis/` found:

1. **Only the HCl proton was quantum, in both isomers.** `quantum_nuc=[0]`
   in the HF-donor inputs and `[2]` in the HCl-donor inputs both select the
   hydrogen bonded to chlorine: the acceptor's proton in F–H···Cl–H, the
   donor proton in Cl–H···F–H. The HF proton was classical in both. Because
   the same molecule's proton was quantum in both isomers, the April energy
   difference is meaningful as a first estimate: the HF-donor form is lower
   by 0.24 to 0.38 kcal/mol for the pure functionals and 0.07 to 0.15 for
   B97, at every basis. These inputs use `QUANTUM_NUC = ['H']`, both
   protons, which adds the missing HF-proton contribution; since that
   proton is the red-shifted donor in the lower isomer, the gap is expected
   to widen slightly.
2. **The chlorine sum-rule defect dominated the raw classification.** At
   grid level 3 the chlorine row was off by 0.07 to 0.13 Eh/Bohr² (aVDZ),
   0.016 to 0.15 (aVTZ) and 0.20 to 0.42 (aVQZ), the same values as in
   HCN···HCl. At DZ and QZ a displaced translation leaked into the lowest
   intermolecular mode as a spurious 54i to 307i cm⁻¹ mode in all 20 runs;
   at TZ it shifted the low modes by up to 70 cm⁻¹ without going imaginary.
   After the diagonal-block repair all 30 spectra have six real modes and
   the repaired lowest mode agrees across bases to a few cm⁻¹ per
   functional. Both isomers are minima everywhere.
3. **Grid level 3 and the loose optimizer settings**, as in the other
   April sets. The starting coordinates were genuinely in Bohr.

## Settings

Those of the other sets (see `../hcn_hf/README.md` and
`../hcn_hbr/README.md`): GAU_TIGHT criteria with `subfrctor=2`, SCF
`conv_tol` 1e-11 / `conv_tol_grad` 1e-6, grid level 5 with grid response
on, 200-step cap, `STOP_IF_OPT_FAILS`, CPHF with 300 Krylov iterations and
a level-shift retry, the sum-rule repair with the raw spectrum printed
alongside.

**Dense halogen radial grids.** `ATOM_GRIDS = {'F': (400, 770), 'Cl':
(400, 770)}` puts 400 radial shells on both halogens, where the level-5
table gives F 105 and Cl 110; hydrogen keeps the table. The HCN···HBr probe
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
`.../dimers/pyscf-master/dimer-alpha-inputs/cneodft/hfhcl`, writing all
logs and results to `.../dimer-alpha-outputs/cneodft/hfhcl`. The tiers are
those of the other sets; with four atoms and at most 247 electronic AOs
they are generous by more than an order of magnitude. Regenerate the
inputs with `python make_inputs.py [date]` after editing the template in
`make_inputs.py`.

## april_analysis/

`analyze_hfcl.py` (with its parser `reanalyze_hessians2.py`) re-diagonalizes
the Hessians printed in the April logs with the six-mode projection, reads
the bonding from the distances, checks the sum rule per atom, applies the
repair and classifies each run on the raw and on the repaired Hessian;
`hf_hcl_classification.txt` is its output for the 30 non-5Z logs. Run it as
`python analyze_hfcl.py <April log dir>`.
