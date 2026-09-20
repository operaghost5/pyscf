# HCN···HBr — conventional DFT re-optimization, Hessian and harmonic frequencies

The all-classical-nuclei counterpart of `../hcn_hbr/`. Forty-five
self-contained PySCF drivers: five functionals (PBE, PW91, BP86, BLYP, B97)
× three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ) × three
isomers of the hydrogen cyanide / hydrogen bromide dimer, conventional
Kohn–Sham DFT (`pyscf.dft.RKS`), no quantum nuclei, bromine all-electron,
400 radial shells on bromine as in `../hcn_hbr/`. Each input

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
python HCNHBr.dft-pbe.global-min.aug-cc-pvdz.2026-09-20.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HCNHBr.dft-<xc>.<global-min|local-min|xlocal-min>.<electronic basis>.<date>.py
```

No `pb4d` token, since there is no protonic basis.

| Token | Isomer | February 2026 DFT result (aVDZ and aVQZ; no aVTZ Hessian survived) |
|---|---|---|
| `global-min` | linear H–C≡N···H–Br, HBr the hydrogen-bond donor | minimum for every functional and basis; N···H 1.99 to 2.18 Å |
| `local-min` | bent N≡C–H···Br–H, the C–H the donor | minimum for every functional and basis, 2.0 to 2.8 kcal/mol up; H···Br 2.81 to 3.01 Å, C–H···Br 174 to 178° |
| `xlocal-min` | linear H–C≡N···Br–H halogen bond, Br the donor | minimum for every functional and basis once the bromine defect is repaired, 2.5 to 3.3 kcal/mol up; N···Br 3.24 to 3.52 Å |

Atom order is H0 C1 N2 H3 Br4 in every file, the same as in `../hcn_hbr/`.

## Starting geometries

`geometries/hcn_hbr_dft_start_geometries.xyz` holds the 45 starting points:
the final geometries of the February 2026 conventional-DFT runs of the same
functional, basis and isomer, with the run's final energy and convergence
status in each comment line. All 45 February optimizations converged (7 to
12 steps for the linear hydrogen-bonded isomer, 13 to 24 for the bent one,
8 to 17 for the halogen-bonded one) and every functional found all three
complexes bound, so no start is substituted. `make_geometries.py`
regenerates the file from the February logs; it reads only the logs whose
method token is `dft`, and the aug-cc-pVTZ geometries come from the
optimizer's own final printout because those logs end before the summary.

## What the February 2026 DFT runs had, and what changed

The classification in `february_analysis/` (`hcn_hbr_dft_classification.txt`)
found:

1. **Every aug-cc-pVTZ Hessian died**, all 15 at the same point: the log
   ends right after the `RHF partial hessian` timing line (750 to 840 s),
   before the exchange-correlation second derivatives and the CPHF, with
   no traceback in the log. The February CNEO-DFT aug-cc-pVTZ runs of this
   complex died in the same tier, at `BEGIN HESSIAN CALCULATION`. Thirty
   runs killed at one stage in one tier points to the TZ job's memory or
   walltime limit, not to the chemistry; the `.err` files or `sacct` for
   those jobs would settle it. The 30 aug-cc-pVDZ and aug-cc-pVQZ Hessians
   finished.
2. **All 30 are minima on the repaired Hessian.** Raw, the BP86 and BLYP
   halogen-bonded aug-cc-pVQZ runs looked like first-order saddle points
   (31.5i and 20.5i cm⁻¹): even with the level-9 grid (200 radial shells)
   the bromine sum-rule violation at aug-cc-pVQZ was 0.011 to 0.025
   Eh/Bohr², enough to displace the translations onto the 15 to 40 cm⁻¹
   modes of these floppy complexes and move them by up to 57 cm⁻¹. At
   aug-cc-pVDZ the violation was only 4e-6 to 4e-4, the same basis
   dependence seen in the CNEO logs (aVTZ worst, aVDZ mildest). The
   repaired lowest modes agree across the two bases.
3. **The pre-2.14 rotor bug hit 23 of the 30 logged frequency lists**,
   dropping one component of a degenerate bend for the linear isomers.
4. **The optimizer settings** were the custom set of all February runs
   (grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3 Å, ΔE 1e-9, 1500-step
   cap) with default SCF tolerances and no grid response.

The energy ordering is unambiguous here because no proton is quantum: the
bent C–H···Br isomer lies 2.4–2.8 (PBE, PW91, BP86) and 2.0–2.2 (BLYP, B97)
kcal/mol above the linear hydrogen-bonded one, and the halogen-bonded
isomer 3.0–3.3 and 2.5–2.7 kcal/mol; both gaps shrink by about 0.2
kcal/mol from aVDZ to aVQZ.

## CNEO-DFT minus DFT (February runs, repaired Hessians, H 1.008 u on both sides)

For the five linear N···H–Br aug-cc-pVDZ pairs, CNEO-DFT lowers the H–Br
stretch by 135 to 149 cm⁻¹ and the C–H stretch by 132 to 134 cm⁻¹, raises
the H–Br libration pair by 35 to 40 cm⁻¹ and the intermolecular stretch by
10 to 12 cm⁻¹, lowers the HCN bend pair by 14 to 17 cm⁻¹, shortens N···H by
0.081 to 0.089 Å and lengthens H–Br by 0.029 to 0.031 Å. For the six
halogen-bonded pairs (aVDZ all functionals, aVQZ B97) the H–Br stretch drops
by 80 to 95 cm⁻¹ and the C–H stretch by 130 to 136, while N···Br lengthens
by 0.028 to 0.045 Å and every intermolecular mode is lower by 3 to 18 cm⁻¹.
As for HCN···HCl, the quantum protons strengthen the hydrogen bond and
weaken the halogen bond. Part 4 of the classification lists every pair.

## Settings

Those of `../hcn_hbr/`, applied to `dft.RKS`: SCF `conv_tol` 1e-11 /
`conv_tol_grad` 1e-6 / 200 cycles, grid level 5 (H 70 × 590, C/N 105 × 770)
with 400 radial shells × 770 angular points on bromine
(`ATOM_GRIDS = {'Br': (400, 770)}`, the level-5 table would give 120), grid
response in the gradient, geomeTRIC GAU_TIGHT with `convergence_energy`
1e-8 and `subfrctor=2`, 200-step cap, `STOP_IF_OPT_FAILS`, `conv_tol_cpscf`
1e-9, 300 Krylov iterations with the level-shift retry, the sum-rule repair
with the raw spectrum printed alongside. `ATOM_GRIDS = {}` would fall back
to the plain level-5 table.

**Masses.** PySCF's isotope-averaged atomic masses (H 1.008, C 12.011,
N 14.007, Br 79.904 u), the standard convention for conventional DFT and
the one the February logs used; the summary prints them. The CNEO mass
convention of the CNEO-DFT inputs (nuclear mass 1.007276 u for the quantum
protons) is not used in these files; the 0.07 % difference in the hydrogen
mass changes a frequency by at most 0.04 %, about 1 cm⁻¹ at 3000 cm⁻¹.

**Why the dense bromine grid and the repair are kept for conventional
DFT.** PySCF's RKS Hessian holds the grid fixed for GGA and hybrid-GGA
functionals exactly as the NEO Hessian does, so the bromine diagonal-block
defect appears here too: 0.011 to 0.025 Eh/Bohr² at aug-cc-pVQZ even with
200 radial shells, and spurious imaginary modes in the raw spectra. The
HCN···HBr probe (`../hcn_hbr/probe/RESULTS.md`) showed that 400 radial
shells remove the defect and that the repair reproduces the converged
spectrum to 0.1 cm⁻¹. Every log prints the per-atom violation before and
after the repair; the bromine number should sit near 1e-6.

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
`.../dimers/pyscf-master/dimer-alpha-inputs/dft/hcnhbr`, writing all logs
and results to `.../dimer-alpha-outputs/dft/hcnhbr`; change the `dft`
directory token in the script if the DFT runs live elsewhere. The tiers
(dz 16 GB / 8 cores, tz 64 GB / 16 cores, qz 128 GB / 32 cores) are those
of the CNEO-DFT sets; the probe measured about 4 GB and 30 minutes on one
core for the CNEO aug-cc-pVTZ Hessian of this complex, so the TZ tier has
ample headroom for whatever killed the February TZ jobs. The sanity check
imports `pyscf.dft` and `pyscf.hessian.rks` instead of `pyscf.neo`.
Regenerate the inputs with `python make_inputs.py [date]` after editing
the template in `make_inputs.py`.

## february_analysis/

`analyze_hcnhx_dft.py` (with its parser `reanalyze_hessians2.py`), the
same generic HCN···HX classifier as in `../hcn_hcl_dft/`;
`hcn_hbr_dft_classification.txt` is its output for the 45 non-5Z DFT logs
(30 with a Hessian), with CNEO-DFT minus DFT shifts for the 11 runs that
have a Hessian in both sets. Run it as
`python analyze_hcnhx_dft.py <DFT log dir> [<CNEO-DFT log dir>]`.
