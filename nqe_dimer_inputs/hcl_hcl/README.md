# (HCl)₂ — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Thirty self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × two starts for the hydrogen chloride dimer, **both**
hydrogen nuclei quantum with the PB4-D protonic basis, no electron-proton
correlation functional, 400 radial shells on chlorine. Each input

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
   Cartesians, Angstrom), the key distances and angles with a C₂h symmetry
   measure, the single-point energy, the sum-rule violation, the number of
   imaginary modes raw and repaired, all 3N frequencies with the
   translations/rotations tagged, and the vibrational frequencies with
   those modes dropped.

```bash
python HClHCl.cneo-dft-pbe.global-min.aug-cc-pvdz.pb4d.2026-09-19.py
```

Next to the log it writes `<stem>.traj.xyz`, `<stem>.opt.xyz` and
`<stem>.hessian.npz`.

## File naming

```
HClHCl.cneo-dft-<xc>.<global-min|c2h-saddle>.<electronic basis>.pb4d.<date>.py
```

| Token | Structure | Expected result |
|---|---|---|
| `global-min` | bent C_s Cl–H···Cl–H, one HCl the donor, the acceptor nearly perpendicular; Cl···Cl 3.68 to 3.88 Å, H···Cl 2.38 to 2.59 Å, Cl–H···Cl 171 to 175°, acceptor tilt 92 to 95° in the April runs | minimum, 6 real modes |
| `c2h-saddle` | C₂h structure with both H–Cl bonds tilted 50° from the Cl···Cl axis in a trans arrangement, Cl···Cl 3.85 Å, the donor–acceptor interchange saddle point | saddle point, 1 imaginary mode, if the optimizer keeps the symmetry |

Atom order is H0 (donor H), Cl1 (acceptor Cl), H2 (acceptor H), Cl3 (donor
Cl) in every file, the same convention as `../hf_hf/`; H0–Cl3 and Cl1–H2
are the covalent bonds, H0···Cl1 the hydrogen bond. Atoms 0 and 2 are the
quantum protons.

## Starting geometries

`geometries/hcl_hcl_start_geometries.xyz` holds the 30 starting points with
provenance in each comment line; `make_geometries.py` regenerates it from
the April logs.

- The 15 `global-min` starts are the final geometries of the April 2026
  runs of the same functional and basis. All 15 converged (9 to 19 steps;
  B97 50, 1191 and 165) and all are minima on the repaired Hessian (see
  `april_analysis/`). In the April logs the donor was H2 (H2–Cl3 donating
  to Cl1), so the two hydrogens are swapped to put the donor at H0; the
  extractor decides the donor from the distances and records the swap.
  For B97/aug-cc-pVQZ the 17 April rerun is used; the 16 April attempt
  died at `BEGIN HESSIAN CALCULATION`.
- The 15 `c2h-saddle` starts are constructed: Cl···Cl 3.85 Å, tilt 50°,
  r(H–Cl) equal to the acceptor H–Cl length of the April geometry of the
  same functional and basis, exact inversion symmetry to floating-point
  precision. An unconstrained minimizer started there has no
  symmetry-breaking gradient component, so it stays on the C₂h surface and
  converges to the saddle point, whose Hessian then has one imaginary
  mode. If numerical noise breaks the symmetry the run slides down to the
  bent minimum; the printout reports the C₂h asymmetry (largest difference
  between the two H–Cl or the two H···Cl distances), so the outcome is
  unambiguous either way.

## What the April 2026 runs had, and what changed

The classification in `april_analysis/` found:

1. **Only the acceptor proton was quantum.** `quantum_nuc=[0]` selects H0,
   but the distances show H2 donating the hydrogen bond in every run. This
   is the reverse of (HF)₂, where the donor was quantum, so the April
   (HF)₂ and (HCl)₂ spectra are not on the same footing; in (HCl)₂ the
   CNEO-lowered stretch (2700 to 2813 cm⁻¹) is the acceptor's and the
   donor stretch (2721 to 2863) is a classical harmonic value. These
   inputs use `QUANTUM_NUC = ['H']`, both protons.
2. **The second start was never defined.** The `lm_start_geom_hydb`
   placeholder holds `DNE` in every field, still with the `F` symbols of the
   (HF)₂ template, so all 20 local-minimum runs died while parsing the
   geometry (two-line logs).
3. **The chlorine sum-rule defect dominated the raw classification.** At
   grid level 3 the chlorine rows were off by 0.07 to 0.13 Eh/Bohr²
   (aVDZ), 0.016 to 0.15 (aVTZ) and 0.20 to 0.42 (aVQZ). At DZ and QZ a
   displaced translation leaked into the lowest intermolecular mode as a
   spurious 206i to 555i cm⁻¹ mode in all ten runs; at TZ it shifted the
   two lowest modes by up to 45 cm⁻¹. After the diagonal-block repair all
   15 spectra have six real modes and the lowest mode agrees across bases
   to about 3 cm⁻¹ per functional (83, 81, 80 cm⁻¹ for PBE).
4. **Grid level 3 and the loose optimizer settings**, as in the other
   April sets. The starting coordinates were genuinely in Bohr.

## Settings

Those of the other sets (see `../hcn_hf/README.md` and
`../hcn_hbr/README.md`): GAU_TIGHT criteria with `subfrctor=2`, SCF
`conv_tol` 1e-11 / `conv_tol_grad` 1e-6, grid level 5 with grid response
on, 200-step cap, `STOP_IF_OPT_FAILS`, CPHF with 300 Krylov iterations and
a level-shift retry, the sum-rule repair with the raw spectrum printed
alongside.

**Dense chlorine radial grid.** `CL_ATOM_GRID = (400, 770)` puts 400 radial
shells on chlorine where the level-5 table gives 110; hydrogen keeps the
table. The HCN···HBr probe (`../hcn_hbr/probe/RESULTS.md`) showed that the
halogen defect is radial quadrature error that 400 shells remove and that
the repair reproduces the converged spectrum to 0.1 cm⁻¹. Every log prints
the per-atom violation before and after the repair; both chlorine numbers
should sit near 1e-6. `CL_ATOM_GRID = None` restores the plain level-5
grid.

## Frequency analysis

All 3N modes with the smallest |f| tagged `TR`, then the same list with
the TR modes dropped, where the number of TR modes is 5 for a linear rotor
and 6 otherwise by PySCF's rotor test at the optimized geometry (6 for
both structures); PySCF's projected analysis as a one-line cross-check.
Both tables use the repaired Hessian when the repair is on. Masses are the
CNEO nuclear masses (`mol.mass`).

## Running on the cluster

```bash
./submit_array.slurm
```

submits three arrays (dz/tz/qz on the `requeue` partition, `brorsenk-lab`
account, every walltime at most 47 h) over the inputs in
`.../dimers/pyscf-master/dimer-alpha-inputs/cneodft/hclhcl`, writing all
logs and results to `.../dimer-alpha-outputs/cneodft/hclhcl`. The tiers are
those of the other sets; with four atoms and at most 242 electronic AOs
they are generous by more than an order of magnitude. To run the minima
alone, leave the `c2h-saddle` files out of the input directory. Regenerate
the inputs with `python make_inputs.py [date]` after editing the template
in `make_inputs.py`.

## april_analysis/

`analyze_hclcl.py` (with its parser `reanalyze_hessians2.py`)
re-diagonalizes the Hessians printed in the April logs with the six-mode
projection, reads the donor from the distances, checks the chlorine sum
rule, applies the repair and classifies each run on the raw and on the
repaired Hessian; `hcl_hcl_classification.txt` is its output for the 30
non-5Z logs. Run it as `python analyze_hclcl.py <April log dir>`.
