# (HBr)₂ — conventional DFT re-optimization, Hessian and harmonic frequencies

The all-classical-nuclei counterpart of `../hbr_hbr/`. Thirty self-contained
PySCF drivers: five functionals (PBE, PW91, BP86, BLYP, B97) × three
electronic bases (aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ) × two starts for
the hydrogen bromide dimer, conventional Kohn–Sham DFT (`pyscf.dft.RKS`),
no quantum nuclei, bromine all-electron, 400 radial shells on bromine as in
`../hbr_hbr/`. Each input

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
   distances and angles with a C₂h symmetry measure, the single-point
   energy, the sum-rule violation, the number of imaginary modes raw and
   repaired, all 3N frequencies with the translations/rotations tagged,
   and the vibrational frequencies with those modes dropped.

```bash
python HBrHBr.dft-pbe.global-min.aug-cc-pvdz.2026-09-20.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HBrHBr.dft-<xc>.<global-min|c2h-saddle>.<electronic basis>.<date>.py
```

No `pb4d` token, since there is no protonic basis.

| Token | Structure | Expected result |
|---|---|---|
| `global-min` | bent C_s Br–H···Br–H, one HBr the donor, the acceptor nearly perpendicular; Br···Br 4.00 to 4.25 Å, H···Br 2.55 to 2.81 Å, Br–H···Br 174 to 178°, acceptor tilt 90 to 93° in the April DFT runs | minimum, 6 real modes |
| `c2h-saddle` | C₂h structure with both H–Br bonds tilted 50° from the Br···Br axis in a trans arrangement, Br···Br 4.15 Å, the donor–acceptor interchange saddle point | saddle point, 1 imaginary mode, if the optimizer keeps the symmetry |

Atom order is H0 (donor H), Br1 (acceptor Br), H2 (acceptor H), Br3 (donor
Br) in every file, the convention of `../hbr_hbr/` and `../hf_hf/`; H0–Br3
and Br1–H2 are the covalent bonds, H0···Br1 the hydrogen bond.

## Starting geometries

`geometries/hbr_hbr_dft_start_geometries.xyz` holds the 30 starting points
with provenance in each comment line; `make_geometries.py` regenerates it
from the April logs (it reads only the logs whose method token is `dft`).

- The 15 `global-min` starts are the final geometries of the April 2026 DFT
  runs of the same functional and basis. All 15 converged (14 to 20 steps)
  and all are minima on the repaired Hessian (see `april_analysis/`). In
  the April logs the donor was H2 (H2–Br3 donating to Br1), as in the CNEO
  April runs, so the extractor swaps the two hydrogens and the two bromines
  to put the donor at H0 and records the swap in the comment line.
- The 15 `c2h-saddle` starts are constructed: Br···Br 4.15 Å, tilt 50°,
  r(H–Br) equal to the acceptor H–Br length of the April DFT geometry of
  the same functional and basis, exact inversion symmetry to floating-point
  precision, the same construction as in `../hbr_hbr/`. An unconstrained
  minimizer started there has no symmetry-breaking gradient component, so
  it stays on the C₂h surface and converges to the saddle point; the
  printout reports the C₂h asymmetry so the outcome is unambiguous either
  way.

## What the April 2026 DFT runs had, and what changed

The classification in `april_analysis/` (`hbr_hbr_dft_classification.txt`)
found:

1. **All 15 global-minimum runs are minima on the repaired Hessian, but
   only 10 on the raw one.** Every aug-cc-pVQZ run showed a spurious 47i to
   77i cm⁻¹ mode and the aug-cc-pVTZ runs had their lowest mode displaced
   upward by up to 57 cm⁻¹ (raw 98, 104 and 92 cm⁻¹ for PBE, BP86 and BLYP
   against 50, 47 and 36 repaired), because even at grid level 9 (200
   radial shells) the bromine sum-rule violation was 5e-3 to 2.7e-2
   Eh/Bohr² at these two bases (4e-6 to 1e-4 at aug-cc-pVDZ) and this
   dimer's lowest mode is only 36 to 53 cm⁻¹. The repair changed the
   vibrations by up to 123 cm⁻¹ (median 49). After it, the lowest mode
   agrees across the three bases to within 2 cm⁻¹ for every functional
   (PBE 50.9, 49.7, 49.3; BLYP 38.3, 36.2, 35.9 cm⁻¹). The Hessians
   finished in 12 to 15 CPHF iterations.
2. **The second start was never defined.** The `lm_start_geom_hydb`
   placeholder holds `DNE` in every field, still with the `F` symbols of the
   (HF)₂ template, so all 15 local-minimum runs died while parsing the
   geometry (two-line logs), as in the CNEO set.
3. **The grid was level 9**, not the level 3 of the CNEO-DFT runs, whose
   Hessians (90 shells on Br) carried 1.4 to 18 Eh/Bohr² and whose raw
   spectra were wrong at every basis.
4. **The optimizer settings** were the custom set of all April runs
   (grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3 Å, ΔE 1e-9, 1500-step
   cap) with default SCF tolerances and no grid response. The starting
   coordinates were genuinely in Bohr.

## CNEO-DFT minus DFT (April runs, repaired Hessians, H 1.008 u on both sides)

The April CNEO-DFT (HBr)₂ runs had only the acceptor proton quantum, the
reverse of (HF)₂. Against the DFT runs, CNEO-DFT lowers the acceptor H–Br
stretch by 81 to 98 cm⁻¹ for PBE, PW91 and BP86 and for BLYP at aVDZ, and
the donor stretch, a classical proton in both, by 6 cm⁻¹. For B97 and for
BLYP at aVTZ/aVQZ the lowered acceptor stretch crosses the donor stretch,
so rank matching splits the shift over the two (17 and 84 cm⁻¹ for
BLYP/aVTZ, 38 and 64 for B97/aVQZ); their sum, 87 to 104 cm⁻¹, is the
quantity to compare. The four intermolecular modes change by only −2 to
+5 cm⁻¹, H···Br shortens by 0.009 to 0.013 Å and the acceptor H–Br
lengthens by 0.023 to 0.027 Å: a quantum acceptor proton hardly touches
the hydrogen bond, as in (HCl)₂ and HF/HCl. The new CNEO inputs in
`../hbr_hbr/` have both protons quantum. Part 4 of the classification lists
every pair.

## Settings

Those of `../hbr_hbr/`, applied to `dft.RKS`: SCF `conv_tol` 1e-11 /
`conv_tol_grad` 1e-6 / 200 cycles, grid level 5 (H 70 × 590) with 400
radial shells × 770 angular points on bromine
(`ATOM_GRIDS = {'Br': (400, 770)}`, the level-5 table would give 120), grid
response in the gradient, geomeTRIC GAU_TIGHT with `convergence_energy`
1e-8 and `subfrctor=2`, 200-step cap, `STOP_IF_OPT_FAILS`, `conv_tol_cpscf`
1e-9, 300 Krylov iterations with the level-shift retry, the sum-rule repair
with the raw spectrum printed alongside, and a dissociation warning at
Br···Br above 7 Å. `ATOM_GRIDS = {}` would fall back to the plain level-5
table.

**Masses.** PySCF's isotope-averaged atomic masses (H 1.008, Br 79.904 u),
the standard convention for conventional DFT and the one the April logs
used; the summary prints them. The CNEO mass convention of the CNEO-DFT
inputs (nuclear mass 1.007276 u for the quantum protons) is not used in
these files; the 0.07 % difference in the hydrogen mass changes a frequency
by at most 0.04 %, about 1 cm⁻¹ at 3000 cm⁻¹.

**Why the dense bromine grid and the repair are kept for conventional
DFT.** PySCF's RKS Hessian holds the grid fixed for GGA and hybrid-GGA
functionals exactly as the NEO Hessian does, so the bromine diagonal-block
defect appears here too, and this floppy dimer is the most sensitive case
of the whole set: at 200 radial shells the raw aug-cc-pVQZ spectra still
showed spurious imaginary modes. The HCN···HBr probe showed that 400 radial
shells remove the defect and that the repair reproduces the converged
spectrum to 0.1 cm⁻¹. Every log prints the per-atom violation before and
after the repair; both bromine numbers should sit near 1e-6.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (6 for
both structures); PySCF's projected analysis as a one-line cross-check.
Both tables use the repaired Hessian when the repair is on.

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/dft/hbrhbr`, writing all logs
and results to `.../dimer-alpha-outputs/dft/hbrhbr`; change the `dft`
directory token in the script if the DFT runs live elsewhere. The tiers are
those of the CNEO-DFT sets and are generous by more than an order of
magnitude for a four-atom conventional-DFT calculation. To run the minima
alone, leave the `c2h-saddle` files out of the input directory. Regenerate
the inputs with `python make_inputs.py [date]` after editing the template
in `make_inputs.py`.

## april_analysis/

`analyze_hxhx_dft.py` (with its parser `reanalyze_hessians2.py`), the same
generic (HX)₂ classifier as in `../hf_hf_dft/`; `hbr_hbr_dft_classification.txt`
is its output for the 30 non-5Z DFT logs (15 with a Hessian), with
CNEO-DFT minus DFT shifts in Part 4 and the quantum proton of each CNEO run
noted. Run it as `python analyze_hxhx_dft.py <DFT log dir> [<CNEO-DFT log dir>]`.
