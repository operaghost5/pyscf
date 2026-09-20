# HF/HCl — conventional DFT re-optimization, Hessian and harmonic frequencies

The all-classical-nuclei counterpart of `../hf_hcl/`. Thirty self-contained
PySCF drivers: five functionals (PBE, PW91, BP86, BLYP, B97) × three
electronic bases (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ) × two isomers of
the hydrogen fluoride / hydrogen chloride dimer, conventional Kohn–Sham DFT
(`pyscf.dft.RKS`), no quantum nuclei, 400 radial shells on both fluorine
and chlorine as in `../hf_hcl/`. Each input

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
   distances and angles with the detected donor, the single-point energy,
   the sum-rule violation, the number of imaginary modes raw and repaired,
   all 3N frequencies with the translations/rotations tagged, and the
   vibrational frequencies with those modes dropped.

```bash
python HFHCl.dft-pbe.hf-donor.aug-cc-pvdz.2026-09-20.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HFHCl.dft-<xc>.<hf-donor|hcl-donor>.<electronic basis>.<date>.py
```

No `pb4d` token, since there is no protonic basis.

| Token | Structure | April 2026 DFT result |
|---|---|---|
| `hf-donor` | F–H···Cl–H, HF the donor, HCl roughly perpendicular to the hydrogen bond; H···Cl 2.21 to 2.35 Å, F–H···Cl 172 to 174°, tilt 92 to 95° | minimum for every functional and basis; the lower isomer by 0.3 to 0.7 kcal/mol |
| `hcl-donor` | Cl–H···F–H, HCl the donor, HF tilted; H···F 1.95 to 2.11 Å, Cl–H···F 169 to 171°, tilt 110 to 118° | minimum for every functional and basis |

The April inputs called these `global_min` and `local_min`; the donor-based
names are those of `../hf_hcl/`, and the extractor assigns them from the
distances. Atom order is H0 F1 H2 Cl3 in every file (H0–F1 is HF, H2–Cl3 is
HCl), the same as in `../hf_hcl/`.

## Starting geometries

`geometries/hf_hcl_dft_start_geometries.xyz` holds the 30 starting points,
the final geometries of the April 2026 DFT runs of the same functional,
basis and isomer, with the run's energy and convergence status in each
comment line. All 30 April optimizations converged (10 to 15 steps) and all
30 Hessians finished. The April `global_min` (HF-donor) inputs had the
hydrogens the other way round, H0 on chlorine and H2 on fluorine, as the
CNEO April inputs did; `make_geometries.py` decides from the distances which
hydrogen belongs to which molecule and swaps the two hydrogens where needed,
so that both isomers share one atom order, and the comment line records the
swap. It reads only the logs whose method token is `dft`.

## What the April 2026 DFT runs had, and what changed

The classification in `april_analysis/` (`hf_hcl_dft_classification.txt`)
found:

1. **All 30 runs are minima, raw and repaired.** The Hessians finished in
   12 to 13 CPHF iterations.
2. **The isomer ordering.** With every nucleus classical the two energies
   are directly comparable: F–H···Cl–H lies below Cl–H···F–H by 0.53–0.67
   (PBE, PW91, BP86), 0.49–0.52 (BLYP) and 0.28–0.36 (B97) kcal/mol, with
   the gap widening slightly from aVDZ to aVTZ/aVQZ. The April CNEO runs,
   which had the HCl proton quantum in both isomers, gave a smaller gap
   (0.24–0.38 for the pure functionals, 0.07–0.15 for B97), so the quantum
   HCl proton favours the isomer in which it donates by about 0.2–0.3
   kcal/mol; the `../hf_hcl/` set with both protons quantum will complete
   the picture.
3. **The grid was level 9** (200 radial shells and 1454 angular points on
   every atom), not the level 3 of the CNEO-DFT runs. The chlorine sum-rule
   violation stayed below 6e-5 Eh/Bohr², the fluorine one below 2.1e-4, and
   the repair moved no vibration by more than 9 cm⁻¹ (median 0.5). The
   April CNEO Hessians at level 3 carried 0.07 to 0.42 Eh/Bohr² on chlorine
   and their raw aVDZ and aVQZ spectra showed a spurious 54i to 307i mode.
4. **The optimizer settings** were the custom set of all April runs
   (grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3 Å, ΔE 1e-9, 1500-step
   cap) with default SCF tolerances and no grid response. The starting
   coordinates were genuinely in Bohr.

## CNEO-DFT minus DFT (April runs, repaired Hessians, H 1.008 u on both sides)

The April CNEO runs had the HCl proton quantum in both isomers. In the
HF-donor isomer that proton is the acceptor's: CNEO-DFT lowers the H–Cl
stretch by 101 to 123 cm⁻¹, the H–F stretch (a classical proton in both) by
7 to 8 cm⁻¹, and changes the intermolecular modes by only −4 to +7 cm⁻¹ and
H···Cl by less than 0.01 Å, so a quantum acceptor proton hardly touches the
hydrogen bond. In the HCl-donor isomer the same proton donates: the H–Cl
stretch drops by 126 to 148 cm⁻¹, the intermolecular modes rise by 2 to
44 cm⁻¹, H–Cl lengthens by 0.025 to 0.029 Å and H···F shortens by 0.053 to
0.066 Å. Part 4 of the classification lists every pair.

## Settings

Those of `../hf_hcl/`, applied to `dft.RKS`: SCF `conv_tol` 1e-11 /
`conv_tol_grad` 1e-6 / 200 cycles, grid level 5 (H 70 × 590) with 400
radial shells × 770 angular points on both halogens
(`ATOM_GRIDS = {'F': (400, 770), 'Cl': (400, 770)}`; the level-5 table would
give F 105 and Cl 110), grid response in the gradient, geomeTRIC GAU_TIGHT
with `convergence_energy` 1e-8 and `subfrctor=2`, 200-step cap,
`STOP_IF_OPT_FAILS`, `conv_tol_cpscf` 1e-9, 300 Krylov iterations with the
level-shift retry, the sum-rule repair with the raw spectrum printed
alongside. Each input checks at the end that H0–F1 and H2–Cl3 are still
the covalent bonds. `ATOM_GRIDS = {}` would fall back to the plain level-5
table.

**Masses.** PySCF's isotope-averaged atomic masses (H 1.008, F 18.998,
Cl 35.45 u), the standard convention for conventional DFT and the one the
April logs used; the summary prints them. The CNEO mass convention of the
CNEO-DFT inputs (nuclear mass 1.007276 u for the quantum protons) is not
used in these files; the 0.07 % difference in the hydrogen mass changes a
frequency by at most 0.04 %, about 1 cm⁻¹ at 3000 cm⁻¹.

**Why the dense halogen grids and the repair are kept for conventional
DFT.** PySCF's RKS Hessian holds the grid fixed for GGA and hybrid-GGA
functionals exactly as the NEO Hessian does, so the chlorine diagonal-block
defect appears here too; at level 3 it produced the spurious imaginary
modes of the CNEO runs. The HCN···HBr probe showed that 400 radial shells
remove it and that the repair reproduces the converged spectrum to
0.1 cm⁻¹. Every log prints the per-atom violation before and after the
repair; both halogen numbers should sit near 1e-6.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (6 for
both isomers); PySCF's projected analysis as a one-line cross-check. Both
tables use the repaired Hessian when the repair is on.

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/dft/hfhcl`, writing all logs
and results to `.../dimer-alpha-outputs/dft/hfhcl`; change the `dft`
directory token in the script if the DFT runs live elsewhere. The tiers are
those of the CNEO-DFT sets and are generous by more than an order of
magnitude for a four-atom conventional-DFT calculation. Regenerate the
inputs with `python make_inputs.py [date]` after editing the template in
`make_inputs.py`.

## april_analysis/

`analyze_hxhy_dft.py` (with its parser `reanalyze_hessians2.py`), the same
generic HX/HY heterodimer classifier as in `../hf_hbr_dft/`;
`hf_hcl_dft_classification.txt` is its output for the 30 non-5Z DFT logs,
with the isomer ranking in the SUMMARY and CNEO-DFT minus DFT shifts in
Part 4. Run it as `python analyze_hxhy_dft.py <DFT log dir> [<CNEO-DFT log dir>]`.
