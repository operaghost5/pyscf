# HCl/HBr — conventional DFT re-optimization, Hessian and harmonic frequencies

The all-classical-nuclei counterpart of `../hcl_hbr/`. Thirty self-contained
PySCF drivers: five functionals (PBE, PW91, BP86, BLYP, B97) × three
electronic bases (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ) × two isomers of
the hydrogen chloride / hydrogen bromide dimer, conventional Kohn–Sham DFT
(`pyscf.dft.RKS`), no quantum nuclei, both halogens all-electron, 400
radial shells on both chlorine and bromine as in `../hcl_hbr/`. Each input

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
python HClHBr.dft-pbe.hcl-donor.aug-cc-pvdz.2026-09-20.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HClHBr.dft-<xc>.<hcl-donor|hbr-donor>.<electronic basis>.<date>.py
```

No `pb4d` token, since there is no protonic basis.

| Token | Structure | April 2026 DFT result |
|---|---|---|
| `hcl-donor` | Cl–H···Br–H, HCl the donor, HBr roughly perpendicular to the hydrogen bond; H···Br 2.51 to 2.74 Å, Cl–H···Br 174 to 177°, tilt 90 to 93° | minimum for every functional and basis; the lower isomer by 0.24 to 0.37 kcal/mol |
| `hbr-donor` | Br–H···Cl–H, HBr the donor, HCl roughly perpendicular; H···Cl 2.44 to 2.69 Å, Br–H···Cl 171 to 176°, tilt 91 to 96° | minimum for every functional and basis once the bromine defect is repaired |

The April inputs called these `global_min` and `local_min`; the donor-based
names are those of `../hcl_hbr/`, and the extractor assigns them from the
distances. Atom order is H0 Cl1 H2 Br3 in every file (H0–Cl1 is HCl, H2–Br3
is HBr), the same as in `../hcl_hbr/`.

## Starting geometries

`geometries/hcl_hbr_dft_start_geometries.xyz` holds the 30 starting points,
the final geometries of the April 2026 DFT runs of the same functional,
basis and isomer, with the run's energy and convergence status in each
comment line. All 30 April optimizations converged (9 to 21 steps) and all
30 Hessians finished. `make_geometries.py` regenerates the file from the
April logs; it reads only the logs whose method token is `dft`, decides
from the distances which molecule donates, and swaps the hydrogens if a
log has them the other way round (none did here).

## What the April 2026 DFT runs had, and what changed

The classification in `april_analysis/` (`hcl_hbr_dft_classification.txt`)
found:

1. **All 30 runs are minima on the repaired Hessian.** The Hessians
   finished in 12 to 14 CPHF iterations. Raw, the BLYP/aug-cc-pVQZ
   HBr-donor run showed a 7i cm⁻¹ mode: even with the level-9 grid (200
   radial shells) the bromine sum-rule violation was 5e-3 to 2.7e-2
   Eh/Bohr² at aug-cc-pVTZ and pVQZ (4e-6 to 1e-4 at pVDZ), and the repair
   moved the soft intermolecular modes of these complexes, 24 to 83 cm⁻¹,
   by up to 52 cm⁻¹ (median 13). The repaired lowest modes are consistent
   across bases where the raw ones scatter. The April CNEO Hessians at
   level 3 carried 1.4 to 18 Eh/Bohr² on bromine and a spurious ~700i mode
   made all nine raw aVQZ spectra look like saddle points.
2. **The isomer ordering is settled.** With every nucleus classical the
   two energies are directly comparable, which the April CNEO energies were
   not (a different proton was quantum in each isomer): Cl–H···Br–H lies
   below Br–H···Cl–H by 0.33–0.37 (PBE, PW91), 0.30–0.31 (BP86), 0.26–0.29
   (BLYP) and 0.24–0.28 (B97) kcal/mol, nearly independent of basis. The
   CNEO re-optimization set in `../hcl_hbr/`, with both protons quantum,
   will show whether the quantum protons change this small gap.
3. **The grid was level 9**, not the level 3 of the CNEO-DFT runs; the
   chlorine violation stayed below 7e-5 Eh/Bohr².
4. **The optimizer settings** were the custom set of all April runs
   (grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3 Å, ΔE 1e-9, 1500-step
   cap) with default SCF tolerances and no grid response. The starting
   coordinates were genuinely in Bohr.

## CNEO-DFT minus DFT (April runs, repaired Hessians, H 1.008 u on both sides)

The April CNEO runs had only each isomer's donor proton quantum. In the
HCl-donor isomer CNEO-DFT lowers the H–Cl stretch by 121 to 150 cm⁻¹,
leaves the H–Br stretch (a classical proton in both) within 4 cm⁻¹, raises
the four intermolecular modes by 5 to 31 cm⁻¹, lengthens H–Cl by 0.025 to
0.030 Å and shortens H···Br by 0.052 to 0.077 Å. In the HBr-donor isomer
it lowers the H–Br stretch by 100 to 124 cm⁻¹, the H–Cl stretch by 2 to
3 cm⁻¹, raises the intermolecular modes by up to 34 cm⁻¹, lengthens H–Br by
0.025 to 0.030 Å and shortens H···Cl by 0.056 to 0.079 Å. In both isomers
the quantum donor proton strengthens the hydrogen bond. Part 4 of the
classification lists every pair; the B97/aug-cc-pVQZ HCl-donor pair uses
the September 2026 CNEO rerun of that case.

## Settings

Those of `../hcl_hbr/`, applied to `dft.RKS`: SCF `conv_tol` 1e-11 /
`conv_tol_grad` 1e-6 / 200 cycles, grid level 5 (H 70 × 590) with 400
radial shells × 770 angular points on both halogens
(`ATOM_GRIDS = {'Cl': (400, 770), 'Br': (400, 770)}`; the level-5 table
would give Cl 110 and Br 120), grid response in the gradient, geomeTRIC
GAU_TIGHT with `convergence_energy` 1e-8 and `subfrctor=2`, 200-step cap,
`STOP_IF_OPT_FAILS`, `conv_tol_cpscf` 1e-9, 300 Krylov iterations with the
level-shift retry, the sum-rule repair with the raw spectrum printed
alongside. Each input checks at the end that H0–Cl1 and H2–Br3 are still
the covalent bonds. `ATOM_GRIDS = {}` would fall back to the plain level-5
table.

**Masses.** PySCF's isotope-averaged atomic masses (H 1.008, Cl 35.45,
Br 79.904 u), the standard convention for conventional DFT and the one the
April logs used; the summary prints them. The CNEO mass convention of the
CNEO-DFT inputs (nuclear mass 1.007276 u for the quantum protons) is not
used in these files; the 0.07 % difference in the hydrogen mass changes a
frequency by at most 0.04 %, about 1 cm⁻¹ at 3000 cm⁻¹.

**Why the dense halogen grids and the repair are kept for conventional
DFT.** PySCF's RKS Hessian holds the grid fixed for GGA and hybrid-GGA
functionals exactly as the NEO Hessian does, so the bromine diagonal-block
defect appears here too: up to 2.7e-2 Eh/Bohr² even with 200 radial
shells, enough to move the softest modes of this floppy complex by tens of
wavenumbers. The HCN···HBr probe showed that 400 radial shells remove it
and that the repair reproduces the converged spectrum to 0.1 cm⁻¹. Every
log prints the per-atom violation before and after the repair; both
halogen numbers should sit near 1e-6.

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
`.../dimers/pyscf-master/dimer-alpha-inputs/dft/hclhbr`, writing all logs
and results to `.../dimer-alpha-outputs/dft/hclhbr`; change the `dft`
directory token in the script if the DFT runs live elsewhere. The tiers are
those of the CNEO-DFT sets and are generous by more than an order of
magnitude for a four-atom conventional-DFT calculation. Regenerate the
inputs with `python make_inputs.py [date]` after editing the template in
`make_inputs.py`.

## april_analysis/

`analyze_hxhy_dft.py` (with its parser `reanalyze_hessians2.py`), the same
generic HX/HY heterodimer classifier as in `../hf_hbr_dft/`;
`hcl_hbr_dft_classification.txt` is its output for the 30 non-5Z DFT logs,
with the isomer ranking in the SUMMARY and CNEO-DFT minus DFT shifts in
Part 4. Run it as `python analyze_hxhy_dft.py <DFT log dir> [<CNEO-DFT log dir>]`.
