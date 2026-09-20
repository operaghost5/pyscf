# HCN···HCl — conventional DFT re-optimization, Hessian and harmonic frequencies

The all-classical-nuclei counterpart of `../hcn_hcl/`. Forty-five
self-contained PySCF drivers: five functionals (PBE, PW91, BP86, BLYP, B97)
× three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ) × three
isomers of the hydrogen cyanide / hydrogen chloride dimer, conventional
Kohn–Sham DFT (`pyscf.dft.RKS`), no quantum nuclei. Each input

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
   distances and angles with the structure type read from the shortest
   intermolecular contact and the off-axis distance, the single-point
   energy, the sum-rule violation, the number of imaginary modes raw and
   repaired, all 3N frequencies with the translations/rotations tagged,
   and the vibrational frequencies with those modes dropped.

```bash
python HCNHCl.dft-pbe.global-min.aug-cc-pvdz.2026-09-20.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HCNHCl.dft-<xc>.<global-min|local-min|xlocal-min>.<electronic basis>.<date>.py
```

No `pb4d` token, since there is no protonic basis.

| Token | Isomer | February 2026 DFT result |
|---|---|---|
| `global-min` | linear H–C≡N···H–Cl, HCl the hydrogen-bond donor | minimum for every functional and basis; N···H 1.92 to 2.09 Å |
| `local-min` | bent N≡C–H···Cl–H, the C–H the donor | minimum for every functional and basis, 2.7 to 3.6 kcal/mol up; H···Cl 2.64 to 2.83 Å, C–H···Cl 172 to 178° |
| `xlocal-min` | linear H–C≡N···Cl–H halogen bond, Cl the donor | minimum with PBE, PW91, BLYP and B97 (N···Cl 3.31 to 3.85 Å, 3.8 to 4.9 kcal/mol up); BP86 found no bound complex |

Atom order is H0 C1 N2 H3 Cl4 in every file, the same as in `../hcn_hcl/`.

## Starting geometries

`geometries/hcn_hcl_dft_start_geometries.xyz` holds the 45 starting points:
the final geometries of the February 2026 conventional-DFT runs of the same
functional, basis and isomer, with the run's final energy and convergence
status in each comment line. All 45 February optimizations converged
(6 to 13 steps for the linear hydrogen-bonded isomer, 17 to 25 for the
bent one, 7 to 13 for the bound halogen-bonded runs) and all 45 Hessians
finished. `make_geometries.py` regenerates the file from the February logs;
it reads only the logs whose method token is `dft`, and when a case has
two logs it keeps the one that finished on the later date (BP86/aug-cc-pVDZ
`xlocal-min` was rerun on 7 July 2026).

**Substituted starts.** The three BP86 `xlocal-min` starts are the PBE
halogen-bonded geometry of the same basis, the convention of `../hcn_hcl/`.
The February BP86 runs from the halogen-bonded start found no bound
complex: at aug-cc-pVTZ and aug-cc-pVQZ the monomers separated to N···Cl
11.8 and 10.6 Å over 133 and 158 steps (the Hessians of two separated
molecules have four and five near-zero imaginary modes, not a saddle point
of a complex), and at aug-cc-pVDZ the run collapsed to the bent C–H···Cl
minimum after 160 steps (the 6 February attempt had separated after 314
steps). These inputs therefore test the BP86 halogen bond afresh with the
tighter settings; if the monomers separate again the input says so
(`> 6 Å` warning) and stops before the Hessian. BLYP, which found no bound
halogen complex in the CNEO-DFT runs, does bind it here, weakly: N···Cl
3.59 to 3.85 Å with intermolecular modes of 12 to 35 cm⁻¹.

## What the February 2026 DFT runs had, and what changed

The classification in `february_analysis/` (`hcn_hcl_dft_classification.txt`)
found:

1. **All 43 bound complexes are minima, raw and repaired.** The two
   apparent saddle points are the separated BP86 pairs above. No CPHF
   failures, unlike the CNEO-DFT set, where every bent-isomer Hessian died
   at the 100-iteration cap.
2. **The grid was level 9**, 200 radial shells and 1454 angular points on
   every atom, not the level 3 of the CNEO-DFT runs. The chlorine sum-rule
   violation was 2e-7 to 5e-4 Eh/Bohr², against 0.02 to 0.4 at level 3 in
   the CNEO logs, and the repair changed no vibration of a bound complex by
   more than 6 cm⁻¹ (B97 halogen-bonded runs, whose softest modes are 11 to
   45 cm⁻¹); the median change is 0.06 cm⁻¹.
3. **The pre-2.14 rotor bug hit 29 of the 45 logged frequency lists.** For
   the two linear isomers PySCF printed nine frequencies where ten are
   expected whenever its rotor test misfired, dropping one component of a
   degenerate bend; the re-analysis restores it.
4. **The optimizer settings** were the custom set of all February runs
   (grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3 Å, ΔE 1e-9, 1500-step
   cap) with default SCF tolerances and no grid response.

The energy ordering is unambiguous here because no proton is quantum: the
bent C–H···Cl isomer lies 3.2–3.6 (PBE, PW91), 3.1–3.4 (BP86), 2.8–3.0
(BLYP) and 2.7–2.9 (B97) kcal/mol above the linear hydrogen-bonded one, and
the halogen-bonded isomer 4.5–4.9 (PBE, PW91), 3.8–4.2 (BLYP) and 3.9–4.2
(B97).

## CNEO-DFT minus DFT (February runs, repaired Hessians, H 1.008 u on both sides)

For the 12 linear N···H–Cl runs with a Hessian in both sets, CNEO-DFT
lowers the H–Cl stretch by 159 to 192 cm⁻¹ and the C–H stretch by 132 to
142 cm⁻¹, raises the H–Cl libration pair by 37 to 50 cm⁻¹ and the
intermolecular stretch by 6 to 15 cm⁻¹, lowers the HCN bend pair by 14 to
23 cm⁻¹, shortens N···H by 0.068 to 0.092 Å and lengthens H–Cl by 0.028 to
0.034 Å: the quantum proton strengthens the hydrogen bond. For the nine
halogen-bonded runs (PBE, PW91, B97) the H–Cl stretch drops by 97 to 117
cm⁻¹ and the C–H stretch by 131 to 139, but N···Cl lengthens by 0.033 to
0.065 Å and every intermolecular mode is lower by 1 to 19 cm⁻¹: the quantum
protons weaken the halogen bond. Part 4 of the classification lists every
pair.

## Settings

Those of `../hcn_hcl/`, applied to `dft.RKS`: SCF `conv_tol` 1e-11 /
`conv_tol_grad` 1e-6 / 200 cycles, grid level 5 (H 70 × 590, C/N 105 × 770,
Cl 110 × 770) with grid response in the gradient, geomeTRIC GAU_TIGHT with
`convergence_energy` 1e-8 and `subfrctor=2`, 200-step cap,
`STOP_IF_OPT_FAILS`, `conv_tol_cpscf` 1e-9, 300 Krylov iterations with the
level-shift retry, the sum-rule repair with the raw spectrum printed
alongside. `ATOM_GRIDS = {}` keeps the level-5 table on every atom to match
`../hcn_hcl/`; `ATOM_GRIDS = {'Cl': (400, 770)}` matches
`../hcn_hcl_densegrid/` and removes the chlorine defect at source.

**Masses.** PySCF's isotope-averaged atomic masses (H 1.008, C 12.011,
N 14.007, Cl 35.45 u), the standard convention for conventional DFT and
the one the February logs used; the summary prints them. The CNEO mass
convention of the CNEO-DFT inputs (nuclear mass 1.007276 u for the quantum
protons) is not used in these files; the 0.07 % difference in the hydrogen
mass changes a frequency by at most 0.04 %, about 1 cm⁻¹ at 3000 cm⁻¹.

**Why the sum-rule repair is kept for conventional DFT.** PySCF's RKS
Hessian holds the grid fixed for GGA and hybrid-GGA functionals exactly as
the NEO Hessian does, so the same chlorine diagonal-block defect appears;
at level 5 it will be far larger than in the level-9 February runs. The
repaired spectra are the ones to use, and the printout shows both.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (5 for the
two linear isomers, 6 for the bent one); PySCF's projected analysis as a
one-line cross-check. Both tables use the repaired Hessian when the repair
is on.

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/dft/hcnhcl`, writing all logs
and results to `.../dimer-alpha-outputs/dft/hcnhcl`; change the `dft`
directory token in the script if the DFT runs live elsewhere. The tiers are
those of the CNEO-DFT sets and are generous for conventional DFT. The
sanity check imports `pyscf.dft` and `pyscf.hessian.rks` instead of
`pyscf.neo`. Regenerate the inputs with `python make_inputs.py [date]`
after editing the template in `make_inputs.py`.

## february_analysis/

`analyze_hcnhx_dft.py` (with its parser `reanalyze_hessians2.py`) serves
every HCN···HX set: it reads the halogen from the logs, re-diagonalizes the
Hessians printed in the February DFT logs with the correct 3N−5 / 3N−6
projection, classifies each structure from its shortest intermolecular
contact, checks the sum rule per atom, applies the repair and classifies
each run on the raw and on the repaired Hessian; given the CNEO-DFT log
directory as a second argument it also prints CNEO-DFT minus DFT frequency
and geometry shifts for the runs with a Hessian in both sets.
`hcn_hcl_dft_classification.txt` is its output for the 46 non-5Z DFT logs
(45 cases). Run it as
`python analyze_hcnhx_dft.py <DFT log dir> [<CNEO-DFT log dir>]`.
