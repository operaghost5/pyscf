# HCN···HF — conventional DFT re-optimization, Hessian and harmonic frequencies

The all-classical-nuclei counterpart of `../hcn_hf/`. Thirty self-contained
PySCF drivers: five functionals (PBE, PW91, BP86, BLYP, B97) × three
electronic bases (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ) × two isomers of
the hydrogen cyanide / hydrogen fluoride dimer, conventional Kohn–Sham DFT
(`pyscf.dft.RKS`), no quantum nuclei. Each input

1. runs a tight DFT SCF at the starting geometry,
2. re-optimizes it with geomeTRIC and the analytic DFT gradient,
3. runs a single point at the optimized geometry,
4. computes the analytic DFT Hessian with the hardened CPHF solve,
5. checks the Hessian's translational sum rule, repairs its diagonal
   blocks, and does the harmonic analysis on both the raw and the
   repaired Hessian, and
6. prints a summary: method, functional, electronic basis, nuclear masses
   and their convention, the DFT grid actually used (level, radial scheme,
   partition, pruning, radial × angular points per element, total after
   pruning), the optimized geometry (symbols and Cartesians, Angstrom), the
   key distances and angles with the detected donor and the off-axis
   distance, the single-point energy, the sum-rule violation, the number
   of imaginary modes raw and repaired, all 3N frequencies with the
   translations/rotations tagged, and the vibrational frequencies with
   those modes dropped.

```bash
python HCNHF.dft-pbe.global-min.aug-cc-pvdz.2026-09-20.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HCNHF.dft-<xc>.<global-min|local-min>.<electronic basis>.<date>.py
```

No `pb4d` token, since there is no protonic basis. `global-min` is the
linear H–C≡N···H–F structure (HF the donor, N the acceptor); `local-min` is
the bent N≡C–H···F–H structure (the C–H the donor, F the acceptor), 4.3 to
5.4 kcal/mol above the linear isomer in the February DFT runs. Atom order
is H0 C1 N2 H3 F4 in every file, the same as in `../hcn_hf/`.

## Starting geometries

`geometries/hcn_hf_dft_start_geometries.xyz` holds the 30 starting points:
the final geometries of the February 2026 conventional-DFT runs of the same
functional, basis and isomer, with the run's final energy and convergence
status in each comment line. All 30 February optimizations converged (6 to
9 steps for the linear isomer, 16 to 29 for the bent one) and all 30
Hessians finished; the CPHF needed 12 to 13 Krylov iterations, against 59
to 69 for the CNEO-DFT Hessians of the same complexes. `make_geometries.py`
regenerates the file from the February logs; it reads only the logs whose
method token is `dft` and ignores the `cneodft` logs that share the
directory.

## What the February 2026 DFT runs had, and what changed

The classification in `february_analysis/` (`hcn_hf_dft_classification.txt`)
found:

1. **All 30 runs are minima, raw and repaired.** The linear complex has ten
   real vibrations (3N−5), the bent complex nine (3N−6). No CPHF failures,
   unlike the CNEO-DFT set, where every bent-isomer and every aug-cc-pVQZ
   Hessian died at the 100-iteration cap.
2. **The grid was level 9**, i.e. 200 radial shells and 1454 angular points
   on every atom, not the level 3 of the CNEO-DFT runs. The fluorine
   sum-rule violation was accordingly small, 4e-8 to 2e-4 Eh/Bohr²
   (hydrogen 1e-8 to 3e-4). Even so, the softest intermolecular mode of the
   bent B97 complex moved by 17 to 26 cm⁻¹ under the repair (raw 85.9, 58.8,
   43.1 cm⁻¹ at DZ/TZ/QZ; repaired 68.2, 75.5, 69.5), and the repaired values
   are the ones that agree across bases. For the other 27 runs the repair
   changed no vibration by more than 6 cm⁻¹, and the median change is 1
   cm⁻¹.
3. **The pre-2.14 rotor bug hit 20 of the 30 logged frequency lists.** For
   the linear complex PySCF printed nine frequencies where ten are expected
   whenever its rotor test misfired, dropping one component of the
   degenerate 78 cm⁻¹ bend; the re-analysis restores it.
4. **The optimizer settings** were the custom set of all February runs
   (grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3 Å, ΔE 1e-9, 1500-step
   cap) with default SCF tolerances and no grid response. Unlike the
   CNEO-DFT runs, none oscillated: the classical gradients at level 9 sit
   below the tight grms threshold.
5. **The halogen-bonded start (`xlocal_min`) was undefined** (`DNE`
   placeholder), so those 15 runs died at setup, as in the CNEO-DFT set. It
   is not part of this set either.

The energy ordering is unambiguous here because no proton is quantum: the
bent isomer lies 5.1–5.4 (PBE, PW91, BP86), 4.8–4.9 (BLYP) and 4.3–4.5
(B97) kcal/mol above the linear one, nearly independent of basis.

## Settings

Those of `../hcn_hf/`, applied to `dft.RKS`: SCF `conv_tol` 1e-11 /
`conv_tol_grad` 1e-6 / 200 cycles, grid level 5 (H 70 × 590, C/N/F 105 ×
770) with grid response in the gradient, geomeTRIC GAU_TIGHT with
`convergence_energy` 1e-8 and `subfrctor=2`, 200-step cap,
`STOP_IF_OPT_FAILS`, `conv_tol_cpscf` 1e-9, 300 Krylov iterations with the
level-shift retry, the sum-rule repair with the raw spectrum printed
alongside. `ATOM_GRIDS = {}` keeps the level-5 table on every atom to match
`../hcn_hf/`; `ATOM_GRIDS = {'F': (400, 770)}` removes the fluorine defect at
source, as in the halogen sets.

**Masses.** The harmonic analysis uses PySCF's isotope-averaged atomic
masses (H 1.008, C 12.011, N 14.007, F 18.998 u), the standard convention
for conventional DFT and the one the February logs used; the summary prints
them. The CNEO mass convention of the CNEO-DFT inputs (nuclear mass
1.007276 u for the quantum protons) is not used in these files. The 0.07 %
difference in the hydrogen mass changes a frequency by at most 0.04 %,
about 1 cm⁻¹ at 3000 cm⁻¹, which is worth stating when CNEO-DFT and DFT
frequencies are compared.

**Why the sum-rule repair is kept for conventional DFT.** PySCF's RKS
Hessian holds the grid fixed for GGA and hybrid-GGA functionals exactly as
the NEO Hessian does (its `grid_response` attribute is not read on that
path), so the same diagonal-block defect appears. At level 5 it will be
larger than in the level-9 February runs; the repaired spectra are the
ones to use, and the printout shows both.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (5 for the
linear isomer, 6 for the bent one); PySCF's projected analysis as a
one-line cross-check. Both tables use the repaired Hessian when the repair
is on.

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/dft/hcnhf`, writing all logs
and results to `.../dimer-alpha-outputs/dft/hcnhf`; change the `dft`
directory token in the script if the DFT runs live elsewhere. The tiers are
those of the CNEO-DFT sets and are generous for conventional DFT: the
February logs show the whole aug-cc-pVQZ run, Hessian included, in well
under an hour single-threaded. The sanity check imports `pyscf.dft` and
`pyscf.hessian.rks` instead of `pyscf.neo`. Regenerate the inputs with
`python make_inputs.py [date]` after editing the template in
`make_inputs.py`.

## february_analysis/

`analyze_hcnhf_dft.py` (with its parser `reanalyze_hessians2.py`)
re-diagonalizes the Hessians printed in the February DFT logs with the
correct 3N−5 / 3N−6 projection, reads the donor from the distances, checks
the sum rule per atom, applies the repair and classifies each run on the
raw and on the repaired Hessian; given the CNEO-DFT log directory as a
second argument it also prints CNEO-DFT minus DFT frequency and geometry
shifts for the ten runs (linear isomer, aVDZ and aVTZ) that have a Hessian
in both sets. `hcn_hf_dft_classification.txt` is its output for the 45
non-5Z DFT logs (30 with a Hessian). Run it as
`python analyze_hcnhf_dft.py <DFT log dir> [<CNEO-DFT log dir>]`.
