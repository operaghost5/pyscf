# HCN···HBr — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Forty-five self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × three isomers of the hydrogen cyanide / hydrogen bromide
dimer, both hydrogen nuclei quantum with the PB4-D protonic basis, no
electron-proton correlation functional, bromine all-electron. Each input

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
python HCNHBr.cneo-dft-pbe.global-min.aug-cc-pvdz.pb4d.2026-09-18.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz` (raw and repaired Hessian, all frequency sets,
geometry, energies).

## File naming

```
HCNHBr.cneo-dft-<xc>.<global-min|local-min|xlocal-min>.<electronic basis>.pb4d.<date>.py
```

| Token | Isomer | February 2026 result |
|---|---|---|
| `global-min` | linear H–C≡N···H–Br, HBr the hydrogen-bond donor | minimum for all functionals, N···H 1.91 to 2.08 Å |
| `local-min` | bent N≡C–H···Br–H, the C–H the donor | stationary point 2.4 to 3.3 kcal/mol up, H(C)···Br 2.75 to 2.94 Å; no aug-cc-pVTZ Hessian survived |
| `xlocal-min` | linear H–C≡N···Br–H halogen bond, Br the donor | bound for every functional, N···Br 3.27 to 3.56 Å, 3.3 to 4.2 kcal/mol up |

Atom order is H, C, N, H, Br in every file; atoms 0 and 3 are the quantum
protons. The energy differences are those of the February runs at their
own final geometries.

## Starting geometries

`geometries/hcn_hbr_start_geometries.xyz` holds the 45 starting points,
the final geometries of the February 2026 CNEO-DFT runs of the same
functional and basis (Angstrom), with the run's final energy and
convergence status in each comment line. Unlike HCN···HCl no substitutions
were needed: every halogen-bonded start is a bound complex. Fourteen of the
February optimizations did not converge, all of them oscillating on the
net-force artefact that `subfrctor=2` removes:

- B97 `global-min`: aug-cc-pVDZ and aug-cc-pVTZ hit the 1500-step cap,
  aug-cc-pVQZ was killed at step 212 (its comment line says so; the
  geometry is the one evaluated at that step).
- `xlocal-min`: PBE, BP86, BLYP and B97 at aug-cc-pVDZ hit the cap; PBE,
  PW91, BP86 and BLYP at aug-cc-pVQZ were killed at steps 202 to 212;
  PW91/aug-cc-pVDZ and PBE and PW91/aug-cc-pVTZ converged only after 900
  to 1450 steps.

`make_geometries.py` regenerates the file from the February logs
(`python make_geometries.py <log dir> geometries/hcn_hbr_start_geometries.xyz`);
where two logs exist for one case (B97/aug-cc-pVQZ `global-min`, 6 and
7 February) it keeps the run that got further.

## Settings

The optimizer, SCF and grid settings are those established for HCN···HF
and HCN···HCl (see `../hcn_hf/README.md`): GAU_TIGHT criteria with
`subfrctor=2`, SCF `conv_tol` 1e-11 / `conv_tol_grad` 1e-6, grid level 5
with grid response on, 200-step cap, `STOP_IF_OPT_FAILS`, CPHF with 300
Krylov iterations and a level-shift retry. One setting is new.

**Dense bromine radial grid.** PySCF's analytic Hessian differentiates the
basis functions but holds the DFT grid points and weights fixed. Under that
approximation the exchange-correlation Hessian violates the translational
sum rule by the quadrature stiffness of each atom's core density on its
own radial grid: isotropic, in that atom's diagonal block, steeply growing
with nuclear charge. On the February Hessians (grid level 3, 90 radial
shells on Br) the bromine row was off by 1.4 to 2.7 Eh/Bohr² against
≤ 1e-3 for every other atom; the translations came out at 690 to 950 cm⁻¹
and the low intermolecular modes were shifted. The probe in `probe/`
measured, for PBE/aug-cc-pVTZ at the February linear geometry:

| Br grid (radial × Lebedev) | Br sum-rule violation, Eh/Bohr² | raw translations, cm⁻¹ |
|---|---|---|
| 120 × 770 (level 5 default) | 7.3e-1 | 438, 438, 491 |
| 120 × 770, no pruning | 7.3e-1 | 438, 438, 491 |
| 240 × 770 | 6.7e-4 | 13, 15, 15 |
| 240 × 770, no pruning | 6.7e-4 | 13, 15, 15 |
| 400 × 974, no pruning | 9.0e-7 | 0.7, 0.7, 0.8 |

Pruning and the angular order do not matter; the radial count does. The
inputs therefore set `BR_ATOM_GRID = (400, 770)` on bromine only (the
other atoms keep the level-5 table), which converges the raw Hessian, and
keep the diagonal-block repair (`SUM_RULE_REPAIR = True`) as a check. The
repair is validated by the same probe: the repaired 120-shell spectrum and
the raw 400-shell spectrum agree to 0.1 cm⁻¹ for every vibration
(43.3, 43.3, 105.9, 441.0, 441.0, 710.3, 710.3, 2129.4, 2223.0, 3217.2 cm⁻¹).
Every input prints the per-atom violation before and after the repair, so
the bromine number in each log should sit near 1e-6; if it does not, the
grid did not take. `BR_ATOM_GRID = None` restores the plain level-5 grid.

## Frequency analysis

As for HCN···HF and HCN···HCl: all 3N modes with the smallest |f| tagged
`TR`, then the same list with the TR modes dropped, where the number of TR
modes is 5 for a linear rotor and 6 otherwise by PySCF's rotor test at the
optimized geometry; PySCF's projected analysis as a one-line cross-check.
Both tables use the repaired Hessian when the repair is on. Masses are the
CNEO nuclear masses (`mol.mass`).

## Cost, from the probe

PBE/aug-cc-pVTZ, 197 electronic AOs, 176 912 grid points at level 5:

| | 1 core (Xeon Gold 6138) | 16 cores |
|---|---|---|
| SCF, 14 cycles | 37 s | 7 s |
| Hessian, all stages | 30 min | 4.9 min |
| of which CPHF (28 Krylov cycles) | 8.8 min | 2.0 min |
| of which electronic second derivatives | 18.7 min | 2.3 min |
| peak RSS | 4.0 GB | 4.1 GB |

The 400-shell bromine grid adds about 60 % to the number of grid points;
on the probe the 240-shell grid cost 10 % more Hessian time and the
400 × 974 unpruned grid 130 % more, so expect something in between. The footprint steps from 2.2 GB during the
optimization to 3.7 GB when the Hessian starts, so a job with a memory
limit between those two numbers dies at `BEGIN HESSIAN CALCULATION` with
no Python traceback, which is what the February aug-cc-pVTZ logs show;
whether that, walltime or preemption killed them is a question for
`sacct`, not for the inputs.

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/cneodft/hcnhbr`, writing all
logs and results to `.../dimer-alpha-outputs/cneodft/hcnhbr`. The tiers
are those of HCN···HCl (16 GB / 8 cores, 64 GB / 16 cores, 128 GB / 32
cores, PYSCF_MAX_MEMORY at half the allocation); against the measured
4 GB and 5 minutes for an aug-cc-pVTZ Hessian on 16 cores they carry a
margin of more than an order of magnitude, which the bent isomers'
longer optimizations and slower CPHF will use part of. Regenerate the
inputs with `python make_inputs.py [date]` after editing the template in
`make_inputs.py`.

## probe/

The stage-by-stage Hessian probes that produced the numbers above, with
the sbatch wrapper and a results summary: `hcnhbr_avtz_hessian_probe.py`
(first version, one core), `hcnhbr_avtz_hessian_probe2.py` (memory
heartbeat, signal tracebacks, split stage E, grid test), `run_probe.sbatch`,
`RESULTS.md`.
