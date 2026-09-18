# HCN···HCl — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Forty-five self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × three isomers of the hydrogen cyanide / hydrogen chloride
dimer, both hydrogen nuclei quantum with the PB4-D protonic basis, no
electron-proton correlation functional. Each input

1. runs a tight CNEO-DFT SCF at the starting geometry,
2. re-optimizes it with geomeTRIC and the analytic CNEO-DFT gradient,
3. runs a single point at the optimized geometry,
4. computes the analytic CNEO-DFT Hessian with the hardened CPHF solve,
5. checks the Hessian's translational sum rule, repairs its diagonal
   blocks, and does the harmonic analysis on both the raw and the
   repaired Hessian, and
6. prints a summary: method, functional, quantum nuclei, EPC, electronic
   and nuclear basis, nuclear masses, the optimized geometry (symbols and
   Cartesians, Angstrom), key intermolecular distances, the single-point
   energy, the sum-rule violation, the number of imaginary modes raw and
   repaired, all 3N frequencies with the translations/rotations tagged,
   and the vibrational frequencies with those modes dropped.

```bash
python HCNHCl.cneo-dft-pbe.global-min.aug-cc-pvdz.pb4d.2026-09-18.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz` (raw and repaired Hessian, all frequency sets,
geometry, energies).

## File naming

```
HCNHCl.cneo-dft-<xc>.<global-min|local-min|xlocal-min>.<electronic basis>.pb4d.<date>.py
```

| Token | Isomer | February 2026 result |
|---|---|---|
| `global-min` | linear H–C≡N···H–Cl, HCl the hydrogen-bond donor | minimum for all functionals, N···H 1.86 to 1.99 Å |
| `local-min` | bent N≡C–H···Cl–H, the C–H the donor | stationary point 3.3 to 4.2 kcal/mol up; no Hessian survived |
| `xlocal-min` | linear H–C≡N···Cl–H halogen bond, Cl the donor | bound with PBE, PW91, B97 (N···Cl 3.36 to 3.60 Å, 5.0 to 5.9 kcal/mol up); unbound with BLYP and BP86 |

Atom order is H, C, N, H, Cl in every file; atoms 0 and 3 are the quantum
protons.

## Starting geometries

`geometries/hcn_hcl_start_geometries.xyz` holds the 45 starting points,
the final geometries of the February 2026 CNEO-DFT runs of the same
functional and basis (Angstrom), with the run's final energy and
convergence status in each comment line. Exceptions, all recorded in the
comment lines and the input docstrings:

- B97/aug-cc-pVQZ `global-min` is the geometry at the last completed
  step of a run killed at step 307 while oscillating on a net-force
  artefact.
- The six BLYP and BP86 `xlocal-min` starts are the PBE halogen-bonded
  geometry of the same basis. The February BLYP and BP86 runs from the
  halogen-bonded start found no bound complex: BP86/aug-cc-pVTZ
  dissociated along the axis to N···Cl 12 Å, the others drifted to 10 to
  13 Å and several Ångström off axis until the step cap or the walltime.
  These inputs therefore test the BLYP and BP86 halogen bond afresh with
  the tighter settings; if the monomers separate again the input says
  so (`N...Cl > 6 Å` warning) and stops before the Hessian.

## Settings

The optimizer, SCF and grid settings are those established for HCN···HF
(see `../hcn_hf/README.md` for the reasoning): GAU_TIGHT criteria with
`subfrctor=2`, SCF `conv_tol` 1e-11 / `conv_tol_grad` 1e-6, grid level 5
with grid response on, 200-step cap, `STOP_IF_OPT_FAILS`. Two Hessian
fixes are specific to what the February HCN···HCl data showed.

**CPHF solver.** 17 of the 45 February Hessians crashed after exactly 100
Krylov iterations, where PySCF's solver raises `RuntimeError`: all 15
bent-isomer Hessians and two linear aug-cc-pVQZ ones. The inputs allow
300 iterations, log the Krylov residuals, and retry once with a 0.2 Eh
level shift on the preconditioner; `conv_tol_cpscf` is 1e-9.

**Translational sum-rule repair.** The 25 February Hessians that did
finish violate the translational sum rule on the chlorine row by 0.07 to
0.13 (aVDZ), 0.016 to 0.15 (aVTZ) and 0.20 to 0.42 Eh/Bohr² (aVQZ),
against ≤ 1e-3 for every other atom and ≤ 1.4e-3 for fluorine in the
HCN···HF Hessians. The error is isotropic and sits in the Cl diagonal
block. In the unprojected spectrum the translations come out at 200 to
560 cm⁻¹ imaginary, and after projection the residual leaks into the
intermolecular bends as spurious imaginary modes of 67i to 352i cm⁻¹ at
aVDZ and aVQZ; at aVTZ the low modes are shifted by tens of cm⁻¹ without
going imaginary. Resetting each diagonal block to minus the sum of its
off-diagonal blocks removes the artefact and makes the spectra agree
across the three bases to a few cm⁻¹ (PBE linear bend 56 / 54 / 52 cm⁻¹
at DZ / TZ / QZ). The likely cause is the missing second-order grid
response in the exchange-correlation Hessian, which grows with nuclear
charge; the finer grid should shrink it, and every input prints the
per-atom violation before and after the repair so this can be checked.
`SUM_RULE_REPAIR = False` turns the repair off; the raw spectrum is
printed in either case.

## Frequency analysis

As for HCN···HF: all 3N modes with the smallest |f| tagged `TR`, then
the same list with the TR modes dropped, where the number of TR modes is
5 for a linear rotor and 6 otherwise by PySCF's rotor test at the
optimized geometry; PySCF's projected analysis as a one-line cross-check.
Both tables use the repaired Hessian when the repair is on. Masses are
the CNEO nuclear masses (`mol.mass`).

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/cneodft/hcnhcl`, writing all
logs and results to `.../dimer-alpha-outputs/cneodft/hcnhcl`. The February
aug-cc-pVQZ HCl runs, single-threaded, took about 260 s per SCF and 400 s
per gradient, a quarter more than HCN···HF, so the tiers are the same.
Regenerate the inputs with `python make_inputs.py [date]` after editing
the template in `make_inputs.py`.
