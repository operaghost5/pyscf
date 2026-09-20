# CNEO-DFT hydrogen-bonded dimers: computational details and the Hessian sum-rule correction

Working notes for the computational-details and methods sections. Part A
lists every system and every parameter of the two generations of
calculations (the original February/April 2026 runs and the re-optimization
inputs generated in September 2026). Part B explains, at the level of a
graduate student meeting the problem for the first time, why the analytic
CNEO-DFT Hessians violate the translational sum rule, what that does to the
harmonic frequencies, how the diagonal-block correction repairs it, why the
correction is exact for this defect, and what the dense radial grids do
instead. Part C gives references and suggested wording. Numbers quoted for
the original runs come from the re-analysis of their logs (the
`april_analysis/` and `*_hessian_reanalysis.txt` files in the repository);
numbers for the probe come from `hcn_hbr/probe/RESULTS.md`.

---

## Part A. Computational details

### A.1 Software

| Item | Original runs (Feb/Apr 2026) | Re-optimization inputs (Sep 2026) |
|---|---|---|
| Electronic-structure code | `theorychemyang/pyscf` fork of PySCF (NEO module), pre-2.14 build, Python 3.10.16, GCC 11.2.0 (from the log headers) | same fork at version 2.14.0 (fork master, Sep 2026), Python 3.12.14, NumPy 2.5.2, SciPy 1.18.1 (the environment of the September probe run; confirm from the production logs) |
| Method | CNEO-DFT, `pyscf.neo.CDFT` | same |
| Geometry optimizer | geomeTRIC through `pyscf.geomopt.geometric_solver` | same |
| Hessian | analytic CNEO-DFT Hessian, `neo.CDFT.Hessian()` | same, with the CPHF settings below |
| Harmonic analysis | `pyscf.hessian.thermo.harmonic_analysis` | same object, called with and without projection, plus the sum-rule check and repair in the driver |
| Functional library | libxc through PySCF | same |

The fork version of the original runs is not printed in their logs; record
the commit hash from the installation used in February/April. The
pre-2.14 build has a rotor-classification defect in
`pyscf.hessian.thermo.rotation_const` (see B.7) that the 2.14.0 build
fixes.

### A.2 Electronic-structure model, common to all systems

- **Method.** Constrained nuclear–electronic orbital DFT (CNEO-DFT):
  quantum protons described by protonic orbitals, with the expectation
  value of each quantum proton's position constrained to the classical
  nuclear position, so that a potential-energy surface, analytic gradients
  and analytic Hessians exist. No electron–proton correlation functional
  (`epc=None`).
- **Functionals.** PBE, PW91, BP86, BLYP, B97. In PySCF/libxc the keyword
  `b97` is the hybrid `HYB_GGA_XC_B97` (Becke 1997, 19.43 % exact
  exchange), not a pure GGA; the other four are pure GGAs.
- **Electronic basis sets.** aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ on all
  atoms; aug-cc-pV5Z was run but excluded from the analysis. All atoms,
  including bromine, are all-electron (nonrelativistic); PySCF's
  `aug-cc-pvXz` for Br is the Wilson–Woon–Peterson–Dunning all-electron
  set, no ECP.
- **Protonic basis.** PB4-D (`nuc_basis='pb4d'`) on every quantum proton.
- **Charge and spin.** Neutral closed-shell singlets throughout.
- **Nuclear masses in the harmonic analysis.** In the re-optimization
  inputs the CNEO masses are used (`mol.mass`): quantum protons carry the
  most-common-isotope atomic mass minus the electron mass, i.e. the proton
  mass 1.007276 u; classical nuclei carry PySCF's isotope-averaged atomic
  masses (C 12.011, N 14.007, F 18.998, Cl 35.45, Br 79.904 u). The
  re-analysis of the original logs used the isotope-averaged atomic masses
  for all atoms (H 1.008 u), which is what PySCF's `harmonic_analysis`
  uses by default; the difference is below 0.1 % in any frequency.

### A.3 Systems and structures

Ten CNEO-DFT input sets exist, plus the conventional-DFT reference set of
A.9; the eight distinct systems and their isomers are
listed with the atom order used in the inputs, the structure, and its
status after the analysis of the original runs. "Minimum" and "saddle"
refer to the repaired Hessian (B.4) unless stated otherwise.

| System | Isomer / start label | Structure | Atom order (quantum protons in bold) | Result of the original runs |
|---|---|---|---|---|
| HCN···HF | `global-min` | linear H–C≡N···H–F, HF the donor | **H** C N **H** F | minimum for every functional at aVDZ/aVTZ (10 Hessians); aVQZ and all `local-min` Hessians died in the CPHF (B.8) |
| | `local-min` | bent N≡C–H···F–H, the C–H the donor, 4.8–6.0 kcal/mol above the global minimum | same | stationary point; no Hessian survived |
| | (`xlocal-min`) | halogen-bonded N···F–H | | never defined in the original inputs (`DNE` placeholder); not part of the HF set |
| HCN···HCl | `global-min` | linear H–C≡N···H–Cl | **H** C N **H** Cl | minimum, N···H 1.86–1.99 Å |
| | `local-min` | bent N≡C–H···Cl–H, 3.3–4.2 kcal/mol up | | stationary point; all 15 Hessians died in the CPHF |
| | `xlocal-min` | linear halogen bond H–C≡N···Cl–H, N···Cl 3.36–3.60 Å, 5.0–5.9 kcal/mol up | | minimum with PBE, PW91, B97; BLYP and BP86 found no bound complex (monomers drifted to 10–13 Å); the new inputs restart those from the PBE geometry |
| HCN···HBr | `global-min` | linear H–C≡N···H–Br, N···H 1.91–2.08 Å | **H** C N **H** Br | minimum (aVDZ; aVTZ Hessians all died, B.8) |
| | `local-min` | bent N≡C–H···Br–H, 2.4–3.3 kcal/mol up, H(C)···Br 2.75–2.94 Å | | stationary point; no Hessian survived |
| | `xlocal-min` | linear halogen bond H–C≡N···Br–H, N···Br 3.27–3.56 Å, 3.3–4.2 kcal/mol up | | minimum for every functional (bound for all five) |
| (HF)₂ | `global-min` | bent C_s F–H···F–H; F···F 2.69–2.75 Å, H···F 1.73–1.81 Å, F–H···F 169–170°, acceptor tilt 109–114° | **H**(donor) F **H** F | minimum, 15/15, raw and repaired |
| | `c2h-saddle` (new) | C₂h trans-tilted interchange saddle, constructed start | | not run in April (`local-min` was undefined) |
| HF/HBr | `hf-donor` (April `global_min`) | F–H···Br–H, HBr perpendicular; H···Br 2.32–2.46 Å, tilt 91–92° | **H** F **H** Br | minimum 15/15 (repaired); raw: aVQZ looked like saddles |
| | `hbr-donor` (April `local_min`) | Br–H···F–H, HF tilted; H···F 1.96–2.14 Å, tilt 108–117° | | minimum 15/15 (repaired) |
| HF/HCl | `hf-donor` (April `global_min`) | F–H···Cl–H, HCl perpendicular; H···Cl 2.21–2.34 Å, tilt 92–94° | **H** F **H** Cl (April inputs had the hydrogens the other way round) | minimum 15/15 (repaired); lower than `hcl-donor` by 0.07–0.38 kcal/mol under the April treatment |
| | `hcl-donor` (April `local_min`) | Cl–H···F–H, HF tilted; H···F 1.90–2.05 Å, tilt 109–116° | | minimum 15/15 (repaired) |
| HCl/HBr | `hcl-donor` (April `global_min`) | Cl–H···Br–H, HBr perpendicular; H···Br 2.46–2.66 Å, tilt 90–93° | **H** Cl **H** Br | minimum 14/14 (repaired); B97/aVQZ optimization killed at step 438 |
| | `hbr-donor` (April `local_min`) | Br–H···Cl–H, HCl perpendicular; H···Cl 2.38–2.61 Å, tilt 91–95° | | minimum 15/15 (repaired) |
| (HCl)₂ | `global-min` | bent C_s Cl–H···Cl–H; Cl···Cl 3.68–3.88 Å, H···Cl 2.38–2.59 Å, Cl–H···Cl 171–175°, tilt 92–95° | **H**(donor) Cl **H** Cl | minimum 15/15 (repaired); raw: DZ and QZ looked like saddles |
| | `c2h-saddle` (new) | C₂h interchange saddle, Cl···Cl 3.85 Å, tilt 50° | | not run in April |
| (HBr)₂ | `global-min` | bent C_s Br–H···Br–H; Br···Br 3.99–4.23 Å, H···Br 2.54–2.80 Å, Br–H···Br 174–177°, tilt 89–92° | **H**(donor) Br **H** Br | minimum 15/15 (repaired); raw spectra wrong at every basis |
| | `c2h-saddle` (new) | C₂h interchange saddle, Br···Br 4.15 Å, tilt 50° | | not run in April |

Two further input sets are variants: `hcn_hcl_densegrid` (HCN···HCl with
400 radial shells on Cl, file token `cl400`, otherwise identical) and the
`c2h-saddle` starts of the three homodimers.

### A.4 Quantum nuclei in the original runs

This differed between systems and is the main reason the original energies
cannot be compared across isomers:

| System (original runs) | `quantum_nuc` | Which proton was quantum |
|---|---|---|
| HCN···HF, HCN···HCl, HCN···HBr | not given (fork default: all H) | both protons |
| (HF)₂ | `[0]` | donor proton only |
| (HCl)₂, (HBr)₂ | `[0]` | acceptor proton only (the donor was atom 2) |
| HF/HCl | `[0]` / `[2]` | the HCl proton in both isomers (acceptor's proton in F–H···Cl–H, donor proton in Cl–H···F–H) |
| HF/HBr, HCl/HBr | `[0]` / `[2]` | each isomer's own donor proton |

The re-optimization inputs treat every hydrogen quantum mechanically
(`quantum_nuc=['H']`) in all ten sets.

### A.5 Settings of the original runs (February and April 2026)

| Setting | Value | Notes |
|---|---|---|
| SCF convergence | PySCF defaults: `conv_tol` 1e-9 Eh, `conv_tol_grad` = √`conv_tol` ≈ 3×10⁻⁵ | |
| DFT grid | PySCF level 3, Treutler–Ahlrichs radial grid, Becke partitioning with Treutler radius adjustment, NWChem pruning | radial shells × Lebedev points: H 50×302, C/N/F 75×302, Cl 80×434, Br 90×434 |
| Grid response in the gradient | off (PySCF default) | |
| Optimizer | geomeTRIC, TRIC coordinates, custom criteria: ΔE 1e-9 Eh, grms 3×10⁻⁶ Eh/Bohr, gmax 4.5×10⁻⁴ Eh/Bohr, drms 1.2×10⁻³ Å, dmax 1.8×10⁻³ Å; `maxiter` 1500 | the grms criterion is 100 times tighter than the gmax criterion and 30 times tighter than Gaussian "tight"; it sits below the gradient noise floor of a level-3 grid without grid response, which is why several optimizations oscillated for hundreds of steps (B97 in particular) or hit the 1500-step cap |
| Starting geometries | Cartesian, `unit='Bohr'` | the (HF)₂ starts were Ångström numbers read as Bohr (compressed by 1.89); all optimizations still reached the correct minimum |
| Hessian | analytic CNEO-DFT Hessian; CPHF by PySCF's Krylov solver, `max_cycle` 100 (fork default), `conv_tol_cpscf` 1e-8 | reaching `max_cycle` raises `RuntimeError`; this killed all bent-isomer HCN···HX Hessians and several aug-cc-pVQZ ones (B.8) |
| Harmonic analysis | `thermo.harmonic_analysis(mol, hessian)` with default projection of translations and rotations | the pre-2.14 rotor test could classify a linear complex as non-linear and drop one bending mode (B.7) |
| Verbosity | 10 | |

### A.6 Settings of the re-optimization inputs (September 2026)

Identical for all ten sets except the dense-grid element and the
system-specific labels.

| Setting | Value |
|---|---|
| SCF | `conv_tol` 1e-11 Eh, `conv_tol_grad` 1e-6, `max_cycle` 200 |
| CPHF for the Hessian | `conv_tol_cpscf` 1e-9; Krylov `max_cycle` 300; level shift 0, one automatic retry with level shift 0.2 Eh if the solver raises |
| DFT grid | level 5: H 70×590; C, N, F 105×770; Cl 110×770; Br 120×770; NWChem pruning, Treutler–Ahlrichs radial grid, Becke partitioning |
| Dense radial grid on the heavy atom(s) | `atom_grid = {X: (400, 770)}` for X = F (hf_hf, hf_hbr, hf_hcl), Cl (hcn_hcl_densegrid, hf_hcl, hcl_hbr, hcl_hcl), Br (hcn_hbr, hf_hbr, hcl_hbr, hbr_hbr); hydrogen and C, N keep the level-5 table |
| Grid response in the gradient | on (`grad.grid_response = True`) |
| Optimizer | geomeTRIC, `convergence_set='GAU_TIGHT'` (grms 1×10⁻⁵, gmax 1.5×10⁻⁵ Eh/Bohr; drms 4×10⁻⁵, dmax 6×10⁻⁵ Å) with `convergence_energy` 1e-8 Eh; `subfrctor=2` (net force and torque projected out of every gradient); 200-step cap; the Hessian is skipped if the optimizer does not converge |
| Starting geometries | Ångström; the final geometry of the original run of the same functional, basis and isomer, or a constructed C₂h structure (homodimers), or the PBE geometry (BLYP/BP86 halogen-bonded HCN···HCl starts) |
| Quantum nuclei | all hydrogens, PB4-D |
| Sum-rule check and repair | per-atom violation printed before and after; diagonal-block repair applied (B.4); raw and repaired spectra both printed |
| Frequency analysis | all 3N eigenvalues of the mass-weighted Hessian; the 5 (linear) or 6 (non-linear) smallest in magnitude tagged as translations/rotations, the rotor type taken from PySCF's rotor test at the optimized geometry; PySCF's projected analysis printed as a cross-check |
| Printed summary | method, functional, quantum nuclei, EPC, electronic and nuclear basis, masses, DFT grid (level, radial scheme, partition, pruning, radial × angular points per element, total after pruning), optimized geometry, key distances and angles (donor detection; C₂h asymmetry for homodimers), single-point energy, sum-rule violation, imaginary-mode counts raw and repaired, all 3N frequencies, vibrational frequencies |

### A.7 Analysis protocol applied to the original logs

For each log the Cartesian Hessian printed by PySCF (Eh/Bohr²), the final
geometry (Å) and the atom symbols were parsed. The Hessian was symmetrized
(the asymmetry was ≤ 10⁻¹² in all cases), mass-weighted and diagonalized in
two ways: (i) without projection, giving all 3N modes, and (ii) after
projecting out the translations and the rotations about the centre of mass
(Eckart vectors, Gram–Schmidt orthonormalized; 5 vectors for a linear
complex, detected geometrically, 6 otherwise), giving the 3N−5 or 3N−6
vibrations. The translational sum-rule violation was evaluated per atom as
max over the 3×3 components of |Σ_B H_AB|. The repaired Hessian of B.4 was
built and diagonalized the same way. A run was classified as a minimum
when all projected vibrational frequencies of the repaired Hessian were
real, and as a saddle point of order n when n of them were imaginary; the
verdict from the raw Hessian is reported alongside for transparency.

### A.8 Sum-rule violations found in the original Hessians (Eh/Bohr²)

| Atom | aug-cc-pVDZ | aug-cc-pVTZ | aug-cc-pVQZ | Effect on the raw spectrum |
|---|---|---|---|---|
| H | ≤ 1×10⁻⁴ | | | none |
| C, N | 1–8×10⁻⁴ | | | none visible |
| F (HCN···HF, (HF)₂, HF/HX) | 1×10⁻³ | 1×10⁻³ | 4–9×10⁻³ | low intermolecular modes of (HF)₂ shifted by 12–22 cm⁻¹ at aVQZ |
| Cl (all Cl systems) | 0.07–0.13 | 0.016–0.15 | 0.20–0.42 | spurious imaginary lowest mode at DZ and QZ (54i–555i cm⁻¹), low modes shifted by up to 70 cm⁻¹ at TZ |
| Br (all Br systems) | 1.4–2.7 | 9.5–18 | 2.7–3.3 | spurious 400i–1045i mode at QZ; spurious 690–2330 cm⁻¹ mode at DZ/TZ; lowest real mode pushed out of the list |

After the repair every Hessian of a bound complex gave real vibrational
frequencies, and for each functional the lowest repaired mode agrees across
the three bases to within a few cm⁻¹, whereas the raw values scatter by
tens to hundreds. The only runs with imaginary modes after repair are the
BLYP and BP86 halogen-bonded HCN···HCl runs whose monomers had drifted
apart to 10–13 Å; those geometries are not stationary points of a bound
complex.

### A.9 Conventional-DFT reference calculations (classical nuclei)

The nuclear quantum effect is measured against conventional Kohn–Sham DFT
(`pyscf.dft.RKS`) with every nucleus classical, run in February 2026 with
the same functionals, bases and starting structures. The HCN···HF set has
been analysed so far (30 non-5Z runs with Hessians; the halogen-bonded
`xlocal_min` start was undefined, as in the CNEO set).

| Setting | Original DFT runs (Feb 2026) | Re-optimization inputs (`hcn_hf_dft/`, Sep 2026) |
|---|---|---|
| Method | `dft.KS(mol, xc=...)`, restricted closed shell | `dft.RKS`, same functionals and bases |
| DFT grid | **level 9**: 200 radial shells × 1454 Lebedev points on every atom (the CNEO runs used level 3) | level 5 (H 70×590, C/N/F 105×770), matching the CNEO re-optimization inputs; `ATOM_GRIDS = {'F': (400, 770)}` available |
| SCF, optimizer, CPHF | as the original CNEO runs (A.5): `conv_tol` 1e-9, custom geomeTRIC criteria, 1500-step cap, RKS Hessian with the default 50-iteration CPHF | as the CNEO re-optimization inputs (A.6) |
| Grid response in the gradient | off | on |
| Masses in the harmonic analysis | PySCF isotope-averaged atomic masses (H 1.008 u) | the same isotope-averaged atomic masses; the CNEO mass convention (nuclear mass 1.007276 u for quantum protons) is used only in the CNEO-DFT inputs. The 0.07 % difference in the hydrogen mass changes a frequency by at most 0.04 %, about 1 cm⁻¹ at 3000 cm⁻¹ |

Results of the original DFT HCN···HF runs (`hcn_hf_dft/february_analysis/`):

- All 30 optimizations converged (6–9 steps linear, 16–29 bent) and all 30
  Hessians finished; the CPHF needed 12–13 Krylov iterations against 59–69
  for the CNEO Hessians of the same complexes, so the CPHF failures of the
  CNEO set (B.8) are specific to the coupled electron–proton equations.
- All 30 are minima on the raw and on the repaired Hessian; the bent isomer
  lies 5.1–5.4 (PBE, PW91, BP86), 4.8–4.9 (BLYP) and 4.3–4.5 (B97) kcal/mol
  above the linear one.
- With 200 radial shells the fluorine sum-rule violation was 4×10⁻⁸ to
  2×10⁻⁴ Eh/Bohr² and the hydrogen violation 10⁻⁸ to 3×10⁻⁴ (larger than in
  the CNEO runs, where the quantum protons have no density cusp on the
  hydrogen grid). Even at this level the softest intermolecular mode of the
  bent B97 complex moved by 17–26 cm⁻¹ under the repair (raw 85.9, 58.8,
  43.1 cm⁻¹ at DZ/TZ/QZ; repaired 68.2, 75.5, 69.5), which shows that the
  soft modes of a floppy complex are sensitive to defects far below the
  level that produces imaginary frequencies. The median change over all
  vibrations was 1 cm⁻¹.
- The pre-2.14 rotor bug (B.7) affected 20 of the 30 logged frequency
  lists: for the linear complex PySCF printed nine frequencies instead of
  ten, dropping one component of the degenerate 78 cm⁻¹ bend.
- PySCF's RKS Hessian holds the grid fixed for GGA and hybrid-GGA
  functionals exactly as the NEO Hessian does (its `grid_response`
  attribute is not read on that path; the grid-response code in
  `hessian/rks.py` serves the meta-GGA and VV10 paths only), so the same
  check and repair apply to the conventional calculations.

CNEO-DFT minus DFT for the ten linear HCN···HF runs with Hessians in both
sets (aVDZ and aVTZ, repaired Hessians, H 1.008 u on both sides): the H–F
stretch is lowered by 231–273 cm⁻¹ and the C–H stretch by 134–144 cm⁻¹,
the C≡N stretch by 6–7 cm⁻¹; the degenerate pair near 650–690 cm⁻¹
(mainly the H–F libration) rises by 10–35 cm⁻¹, the pair near 700–770 cm⁻¹
(mainly the HCN bend) falls by up to 24 cm⁻¹, and the intermolecular
stretch rises by 13–20 cm⁻¹; the N···H distance shortens by
0.040–0.051 Å and the H–F bond lengthens by 0.026–0.028 Å. These are the
harmonic CNEO frequencies against harmonic DFT frequencies, i.e. the
quantum-proton correction to the harmonic picture, not an anharmonic
correction to the DFT value.

Results of the original DFT HCN···HCl runs (`hcn_hcl_dft/february_analysis/`,
45 cases, one rerun):

- All 45 optimizations converged and all 45 Hessians finished. The 43
  bound complexes are minima raw and repaired: 15 linear N···H–Cl, 15 bent
  C–H···Cl and 13 linear halogen-bonded N···Cl–H (PBE, PW91, BLYP, B97;
  BLYP binds it weakly, N···Cl 3.59–3.85 Å with intermolecular modes of
  12–35 cm⁻¹). BP86 found no bound halogen complex: the monomers separated
  at aVTZ and aVQZ and the aVDZ run collapsed to the bent minimum, as in
  the CNEO runs, where BLYP failed too. The bent isomer lies 2.7–3.6
  kcal/mol and the halogen-bonded isomer 3.8–4.9 kcal/mol above the linear
  hydrogen-bonded one.
- With 200 radial shells the chlorine violation was 2×10⁻⁷ to 5×10⁻⁴
  Eh/Bohr² (0.02–0.4 at level 3 in the CNEO logs) and the repair changed no
  vibration of a bound complex by more than 6 cm⁻¹ (median 0.06). The rotor
  bug affected 29 of the 45 logged lists.
- CNEO-DFT minus DFT for the 12 linear hydrogen-bonded pairs: H–Cl stretch
  −159 to −192 cm⁻¹, C–H stretch −132 to −142, H–Cl libration pair +37 to
  +50, HCN bend pair −14 to −23, intermolecular stretch +6 to +15; N···H
  shorter by 0.068–0.092 Å, H–Cl longer by 0.028–0.034 Å. For the nine
  halogen-bonded pairs (PBE, PW91, B97): H–Cl stretch −97 to −117, C–H
  −131 to −139, every intermolecular mode lower by 1–19 cm⁻¹ and N···Cl
  longer by 0.033–0.065 Å. The quantum protons strengthen the hydrogen bond
  and weaken the halogen bond.

Results of the original DFT HCN···HBr runs (`hcn_hbr_dft/february_analysis/`,
45 runs, bromine all-electron):

- All 45 optimizations converged and every functional bound all three
  isomers (the halogen-bonded N···Br–H complex included, N···Br 3.24–3.52
  Å). All 15 aug-cc-pVTZ Hessians died at the same point, right after the
  RHF partial Hessian and before the XC second derivatives, with no error
  in the log; the CNEO aug-cc-pVTZ runs of this complex died in the same
  tier at the start of the Hessian, so a TZ-tier memory or walltime limit
  is the likely cause. The 30 aug-cc-pVDZ and aug-cc-pVQZ Hessians
  finished.
- All 30 are minima on the repaired Hessian. Raw, the BP86 and BLYP
  halogen-bonded aVQZ runs looked like first-order saddle points (31.5i and
  20.5i cm⁻¹): even with 200 radial shells the bromine violation at aVQZ
  was 0.011–0.025 Eh/Bohr² and moved the 15–40 cm⁻¹ modes of these floppy
  complexes by up to 57 cm⁻¹, while at aVDZ it was only 4×10⁻⁶ to 4×10⁻⁴.
  The bent isomer lies 2.0–2.8 kcal/mol and the halogen-bonded isomer
  2.5–3.3 kcal/mol above the linear hydrogen-bonded one. The rotor bug
  affected 23 of the 30 logged lists.
- CNEO-DFT minus DFT for the five linear hydrogen-bonded aVDZ pairs: H–Br
  stretch −135 to −149 cm⁻¹, C–H stretch −132 to −134, H–Br libration pair
  +35 to +40, HCN bend pair −14 to −17, intermolecular stretch +10 to +12;
  N···H shorter by 0.081–0.089 Å, H–Br longer by 0.029–0.031 Å. For the six
  halogen-bonded pairs: H–Br stretch −80 to −95, C–H −130 to −136, every
  intermolecular mode lower by 3–18 cm⁻¹ and N···Br longer by 0.028–0.045
  Å. The same pattern as HCN···HCl.

Results of the original DFT (HF)₂ runs (`hf_hf_dft/april_analysis/`, April
2026, 30 runs):

- The 15 global-minimum runs converged (30–42 steps, from the same
  Å-read-as-Bohr compressed start as the CNEO runs) and all 15 Hessians
  finished in 11 CPHF iterations; all 15 are minima raw and repaired. The
  15 local-minimum runs died on the undefined `DNE` geometry, as in the
  CNEO set. Bent C_s structure: F···F 2.71–2.77 Å, H···F 1.77–1.86 Å,
  F–H···F 169–171°, acceptor tilt 109–114°.
- With 200 radial shells the fluorine violation was 1×10⁻⁶ to 2×10⁻⁴
  Eh/Bohr² and the hydrogen violation up to 2.4×10⁻⁴; the repair moved no
  vibration by more than 6 cm⁻¹ (median 0.7).
- CNEO-DFT minus DFT (the April CNEO runs had only the donor proton
  quantum): donor H–F stretch −206 to −240 cm⁻¹, acceptor H–F stretch (a
  classical proton in both) −6 to −7, the four intermolecular modes +6 to
  +38; donor H–F longer by 0.023–0.025 Å, H···F shorter by 0.034–0.046 Å,
  F···F shorter by 0.010–0.022 Å.

Results of the original DFT HF/HBr runs (`hf_hbr_dft/april_analysis/`, April
2026, 30 runs):

- All 30 optimizations converged (7–17 steps) and all 30 Hessians finished
  (12–14 CPHF iterations); all 30 are minima raw and repaired. With every
  nucleus classical the isomer energies are directly comparable, which the
  April CNEO energies were not: F–H···Br–H lies below Br–H···F–H by
  0.96–1.10 (PBE, PW91, BP86), 0.86–0.90 (BLYP) and 0.62–0.69 (B97)
  kcal/mol. Geometries: HF-donor H···Br 2.36–2.50 Å, tilt 91–92°;
  HBr-donor H···F 2.02–2.20 Å, tilt 110–119°.
- The bromine violation was 4×10⁻⁶ to 10⁻⁴ Eh/Bohr² at aVDZ but 5×10⁻³ to
  2.7×10⁻² at aVTZ and aVQZ with 200 radial shells; the repair moved the
  intermolecular modes by up to 12 cm⁻¹ (median 6) without producing any
  spurious imaginary mode, the softest mode (63 cm⁻¹) being stiffer than in
  HCN···HBr.
- CNEO-DFT minus DFT (each isomer's donor proton quantum in the April CNEO
  runs): HF-donor isomer, H–F stretch −203 to −234 cm⁻¹, H–Br stretch (a
  classical proton in both) within 1.5 cm⁻¹, intermolecular modes +3 to
  +26, H–F longer by 0.023–0.025 Å, H···Br shorter by 0.035–0.050 Å;
  HBr-donor isomer, H–Br stretch −102 to −122, H–F stretch −5 to −7,
  intermolecular modes +4 to +45, H–Br longer by 0.025–0.030 Å, H···F
  shorter by 0.058–0.071 Å.

Results of the original DFT HF/HCl runs (`hf_hcl_dft/april_analysis/`, April
2026, 30 runs):

- All 30 optimizations converged (10–15 steps) and all 30 Hessians finished
  (12–13 CPHF iterations); all 30 are minima raw and repaired. The isomer
  energies are directly comparable: F–H···Cl–H lies below Cl–H···F–H by
  0.53–0.67 (PBE, PW91, BP86), 0.49–0.52 (BLYP) and 0.28–0.36 (B97)
  kcal/mol. The April CNEO runs, with the HCl proton quantum in both
  isomers, gave 0.24–0.38 and 0.07–0.15, so the quantum HCl proton favours
  the isomer in which it donates by about 0.2–0.3 kcal/mol. Geometries:
  HF-donor H···Cl 2.21–2.35 Å, tilt 92–95°; HCl-donor H···F 1.95–2.11 Å,
  tilt 110–118°. The April HF-donor inputs had the hydrogens in the other
  order; the extractor reorders them from the distances.
- With 200 radial shells the chlorine violation stayed below 6×10⁻⁵
  Eh/Bohr² and the fluorine violation below 2.1×10⁻⁴; the repair moved no
  vibration by more than 9 cm⁻¹ (median 0.5).
- CNEO-DFT minus DFT (HCl proton quantum in both April CNEO isomers): in
  the HF-donor isomer, where that proton is the acceptor's, the H–Cl
  stretch drops 101–123 cm⁻¹ while the intermolecular modes change by only
  −4 to +7 cm⁻¹ and H···Cl by less than 0.01 Å; in the HCl-donor isomer,
  where it donates, the H–Cl stretch drops 126–148 cm⁻¹, the intermolecular
  modes rise 2–44 cm⁻¹, H–Cl lengthens by 0.025–0.029 Å and H···F shortens
  by 0.053–0.066 Å. A quantum acceptor proton hardly affects the hydrogen
  bond; a quantum donor proton strengthens it.

Results of the original DFT HCl/HBr runs (`hcl_hbr_dft/april_analysis/`,
April 2026, 30 runs):

- All 30 optimizations converged (9–21 steps) and all 30 Hessians finished
  (12–14 CPHF iterations). All 30 are minima on the repaired Hessian; raw,
  the BLYP/aVQZ HBr-donor run showed a 7i mode. The isomer energies are
  directly comparable: Cl–H···Br–H lies below Br–H···Cl–H by 0.33–0.37
  (PBE, PW91), 0.30–0.31 (BP86), 0.26–0.29 (BLYP) and 0.24–0.28 (B97)
  kcal/mol. Geometries: HCl-donor H···Br 2.51–2.74 Å, tilt 90–93°;
  HBr-donor H···Cl 2.44–2.69 Å, tilt 91–96°.
- The bromine violation was 4×10⁻⁶ to 10⁻⁴ Eh/Bohr² at aVDZ but 5×10⁻³ to
  2.7×10⁻² at aVTZ and aVQZ with 200 radial shells; the repair moved the
  24–83 cm⁻¹ intermolecular modes by up to 52 cm⁻¹ (median 13). The
  chlorine violation stayed below 7×10⁻⁵.
- CNEO-DFT minus DFT (each isomer's donor proton quantum in the April CNEO
  runs): HCl-donor isomer, H–Cl stretch −121 to −150 cm⁻¹, H–Br stretch
  within 4 cm⁻¹, intermolecular modes +5 to +31, H–Cl longer by
  0.025–0.030 Å, H···Br shorter by 0.052–0.077 Å; HBr-donor isomer, H–Br
  stretch −100 to −124, H–Cl stretch −2 to −3, intermolecular modes up to
  +34, H–Br longer by 0.025–0.030 Å, H···Cl shorter by 0.056–0.079 Å.

Results of the original DFT (HCl)₂ runs (`hcl_hcl_dft/april_analysis/`,
April 2026, 30 runs):

- The 15 global-minimum runs converged (13–16 steps) and all 15 Hessians
  finished (11–13 CPHF iterations); all 15 are minima raw and repaired. The
  15 local-minimum runs died on the undefined `DNE` geometry, as in the CNEO
  set. Bent C_s structure: Cl···Cl 3.69–3.89 Å, H···Cl 2.39–2.60 Å,
  Cl–H···Cl 172–176°, acceptor tilt 92–95°; the donor was H2 in the April
  logs and the extractor moves it to H0.
- With 200 radial shells the chlorine violation stayed below 10⁻⁴ Eh/Bohr²
  and the repair moved no vibration by more than 0.6 cm⁻¹ (the April CNEO
  Hessians at level 3 carried 0.07–0.42 and showed spurious 206i–555i
  modes).
- CNEO-DFT minus DFT (the April CNEO runs had only the acceptor proton
  quantum): acceptor H–Cl stretch −103 to −119 cm⁻¹ where the two stretches
  stay separated (PBE, PW91, BP86); for BLYP and B97 the lowered acceptor
  stretch crosses the donor stretch and rank matching splits a total shift
  of 108–125 cm⁻¹ over both. The donor stretch (classical in both) moves by
  −8 to −10, the intermolecular modes by −2 to +8; H···Cl shortens by only
  0.010–0.015 Å. A quantum acceptor proton hardly touches the hydrogen
  bond, as in HF/HCl.

---

## Part B. Why the Hessian needed a sum-rule correction, and how it works

### B.1 What the translational sum rule says

Write the Cartesian Hessian in atom blocks, H_AB, each a 3×3 matrix of
second derivatives ∂²E/∂R_A∂R_B. If every atom is displaced by the same
vector **t**, the energy of an isolated molecule cannot change, because
nothing in the Hamiltonian depends on where the molecule sits in space.
Expanding E to second order in that displacement gives

  Σ_A Σ_B **t**ᵀ H_AB **t** = 0 for every **t**,

and applying the same argument to the gradient (∂E/∂R_A does not change
under translation either) gives the stronger, row-wise statement

  Σ_B H_AB = 0 for every atom A (each of the nine components).

This is the translational sum rule. In solid-state language it is the
acoustic sum rule: the uniform translation is an acoustic phonon of zero
frequency. Its consequences for the vibrational analysis are that the
three translation vectors are exact null vectors of the mass-weighted
Hessian, so three of the 3N eigenvalues are exactly zero and the remaining
3N−3 are unaffected by them. Rotations give an analogous condition: at a
stationary point (zero gradient) the three infinitesimal rotations about
the centre of mass are also null vectors, which is why 3N−6 vibrations
remain (3N−5 for a linear molecule, which has only two rotations). Wilson,
Decius and Cross treat this in the language of the Eckart conditions.

A Hessian that satisfies these conditions can be analysed in two
equivalent ways: diagonalize all 3N modes and discard the six (or five)
zero eigenvalues, or project the translation and rotation vectors out first
and diagonalize the 3N−6 dimensional remainder. PySCF's
`harmonic_analysis` does the latter. The two routes agree only when the
Hessian actually obeys the sum rules. When it does not, the projection
still produces 3N−6 numbers, but they are no longer the vibrational
frequencies of the system, because the part of the Hessian that lives in
the translation subspace is not removed by the projection; it is folded
into the remaining modes. This is the mechanism behind every spurious
imaginary frequency in the original data.

### B.2 Why the CNEO-DFT Hessian violates it: fixed-grid quadrature

The exchange–correlation energy in a DFT code is evaluated by numerical
quadrature on a molecular grid built from atom-centred spherical grids
(radial shells × Lebedev angular points) stitched together with Becke's
partition functions. The grid points and weights move with the atoms; when
the molecule translates, the whole grid translates with it and the
quadrature value is exactly invariant.

Analytic derivatives of E_xc therefore contain two kinds of terms: the
derivatives of the density (through the basis functions, which move with
the nuclei) and the derivatives of the grid points and weights. The second
kind is usually called the grid-weight derivative or grid response.
PySCF's gradient can include it (`grid_response=True`); PySCF's analytic
Hessian, including the NEO version, does not include it at all. It
differentiates the basis functions twice while holding every grid point
and weight fixed.

That approximation breaks translational invariance in a specific way.
Consider the XC energy of a single atom's core density integrated on that
atom's own radial grid. If the density is shifted by a small vector **t**
while the grid stays put, the quadrature error changes; to second order it
changes as ½ **t**ᵀ **K**_A **t**, where **K**_A is the curvature of the
quadrature error with respect to sliding the density off the grid. The
fixed-grid Hessian contains exactly this term in the diagonal block H_AA,
and nothing compensates it in the off-diagonal blocks, so the row sum
Σ_B H_AB is no longer zero but equals **K**_A. The quantity we measured as
the "sum-rule violation" of atom A is the largest component of **K**_A.

Three properties of **K**_A follow from this picture and were all
observed:

1. It is isotropic (the radial grid has no preferred direction), so the
   violation is the same in x, y and z.
2. It lives in the diagonal block of that atom only. The off-diagonal
   blocks, which describe how a force on A changes when B moves, are exact
   to the usual precision.
3. It grows steeply with nuclear charge, because the core density becomes
   sharper as Z increases and a fixed set of radial shells resolves it less
   well. Measured values on the level-3 grid: H ≈ 10⁻⁷ to 10⁻⁴, C and N
   ≈ 10⁻⁴, F ≈ 10⁻³, Cl 0.02–0.4, Br 1–18 Eh/Bohr².

The decisive test is the dependence on the radial grid. For HCN···HBr at
the PBE/aug-cc-pVTZ geometry, with everything else fixed, the bromine
violation was 0.73 Eh/Bohr² with 120 radial shells (PySCF level 5),
6.7×10⁻⁴ with 240 shells and 9×10⁻⁷ with 400 shells, roughly three
orders of magnitude at each step, while switching the pruning off or raising the
angular order from 770 to 974 points changed nothing at the 10⁻⁴ level.
A missing physical term would not vanish with a denser grid; a quadrature
artefact does. The SCF energy itself changed by only 0.03 µEh across these
grids, which shows how much more sensitive the fixed-grid second
derivative is than the energy.

Why the basis matters as well: the aug-cc-pVTZ violation for bromine
(9.5–18) is larger than the aug-cc-pVDZ (1.4–2.7) and aug-cc-pVQZ (2.7–3.3)
ones at the same grid. The curvature **K**_A depends on how sharply the
core density is represented, which differs between the contraction schemes
of the three basis sets; the effect is not monotonic in the cardinal
number.

### B.3 What the defect does to the frequencies

An isotropic error δ **1** added to the diagonal block of atom A changes
the mass-weighted Hessian in the translation direction. The three
translation vectors, normalized in mass-weighted coordinates, have
components √(m_B/M) on atom B, where M is the total mass; the Rayleigh
quotient of the defective Hessian along such a vector is δ/M. If nothing
else mixed in, the translations would appear as modes of frequency

  ω_T ≈ √(δ/M).

For HCN···HBr, δ = 0.73 Eh/Bohr² and M = 108 u give ω_T ≈ 423 cm⁻¹; the
raw Hessian indeed showed its three translations at 438, 438 and 491 cm⁻¹
instead of zero. A positive δ produces spurious real modes, a negative δ
spurious imaginary ones; both signs occurred. The chlorine defect at
aug-cc-pVDZ and aug-cc-pVQZ and the bromine defect at aug-cc-pVQZ
displaced the translations to imaginary frequency, the bromine defect at
aug-cc-pVDZ and aug-cc-pVTZ to real frequency. The sign is a property of
the quadrature error of the particular basis and grid, not of the element.

The translations do not stay pure, however. Once they carry a finite
frequency they mix with the genuine low-frequency intermolecular modes
(the intermolecular stretch and the librations, 40–200 cm⁻¹ in these
complexes) whenever the translation frequency is comparable, and after
PySCF's projection the residual is redistributed over the 3N−6
"vibrations". The symptoms depend on the size of δ:

- Small δ (fluorine, 10⁻³): translations at a few cm⁻¹, low intermolecular
  modes shifted by a few cm⁻¹; at 10⁻² (fluorine, aVQZ) the two lowest
  modes of (HF)₂ shift by 12–22 cm⁻¹.
- Intermediate δ (chlorine, 0.1–0.4): the displaced translation lands on
  top of the lowest intermolecular mode, and the projected spectrum acquires
  one imaginary frequency of 50–550i cm⁻¹ that is not a real curvature.
  This is why every aug-cc-pVDZ and aug-cc-pVQZ Hessian of the chlorine
  systems looked like a first-order saddle point.
- Large δ (bromine, 1–18): the displaced translation appears as a spurious
  mode anywhere from 400 to 2300 cm⁻¹, the lowest real mode is pushed out
  of the 3N−6 list altogether, and at aug-cc-pVQZ a spurious 700–1045i
  mode appears. For (HBr)₂ the raw vibrational list was wrong at every
  basis, not only where an imaginary mode showed.

None of this affects the geometry optimizations, which use gradients, or
the high-frequency intramolecular stretches, which are far from the
translation frequency and change by less than 1 cm⁻¹ under the repair.

### B.4 The repair: enforcing the sum rule on the diagonal blocks

Because the defect sits in the diagonal blocks and the off-diagonal blocks
are correct, the sum rule determines the correct diagonal blocks uniquely:

  H_AA ← − Σ_{B≠A} H_AB, followed by one symmetrization H ← ½(H + Hᵀ).

This replaces each atom's self-interaction block by minus the sum of its
interactions with every other atom, which is what translational invariance
demands, and leaves everything else untouched. Applied to the HCN···HBr
probe Hessian computed on the 120-shell grid, the repaired spectrum
reproduced the spectrum of the raw 400-shell Hessian to 0.1 cm⁻¹ in every
one of the ten vibrations (43.3, 43.3, 105.9, 441.0, 441.0, 710.3, 710.3,
2129.4, 2223.0, 3217.2 cm⁻¹). That agreement is the validation: the repair
recovers what a converged grid gives, without the converged grid.

Two consistency checks come with it. First, the repaired translations sit
at 0 to 1 cm⁻¹ and the rotations at 0 to 3 cm⁻¹, as they should. Second,
for a given functional the repaired lowest intermolecular mode agrees
across aug-cc-pVDZ, TZ and QZ to a few cm⁻¹ (for example 83, 81 and 80
cm⁻¹ for PBE (HCl)₂; 52, 51 and 50 cm⁻¹ for PBE (HBr)₂), whereas the raw
values are scattered or imaginary. Basis-set convergence of a low
intermolecular mode is a strong test because the defect depends on the
basis and the grid in an uncorrelated way.

This correction is the standard "simple" acoustic-sum-rule enforcement of
lattice-dynamics codes (for instance the `asr='simple'` option of Quantum
ESPRESSO's `matdyn`/`dynmat`), where the self force constant of each atom
is reset so that the force constants sum to zero. Gonze and Lee discuss the
same violation arising from real-space grids in DFPT phonon calculations
and the same fix.

What the repair cannot do: it corrects only a diagonal-block error. An
error distributed over off-diagonal blocks, or one that depended on
direction, would need the rotational sum rule and a least-squares
projection instead. We tested for such errors indirectly: the rotational
modes came out at 0 to 3 cm⁻¹ after the diagonal repair alone, and the
cross-basis agreement held, so nothing beyond the diagonal blocks was
detectably wrong.

### B.5 Why projecting the translations out does not work instead

PySCF's `harmonic_analysis` projects the six translation and rotation
vectors out of the mass-weighted Hessian before diagonalizing. It is
tempting to think this removes the defect. It does not, for a simple
reason: the projection assumes the Hessian is already exactly invariant,
so that the translation subspace contains nothing but numerical noise. If
the Hessian has a component δ/M along the translation direction, the
projection removes that component's action within the translation subspace
but not its coupling to the vibrational subspace; the coupling terms are
what mix the displaced translation into the low modes. Mathematically,
P H P with P the projector onto the vibrational subspace is not the
vibrational block of the exact Hessian when H has translation–vibration
cross terms, and those cross terms are generated by the same diagonal-block
error. Repairing the diagonal blocks first makes P H P correct; projecting
alone leaves the contamination in place. This is why the original
PySCF-projected frequencies were wrong even though PySCF had "removed" the
translations.

### B.6 The real fix: denser radial grids

The repair is a correction; the dense radial grid is the cure. Because the
defect is quadrature stiffness of the core density on the atom's own
radial grid, more radial shells on that atom remove it at the source. The
re-optimization inputs place 400 radial shells on each halogen (F, Cl, Br)
with 770 angular points, against 105, 110 and 120 shells in PySCF's level-5
table, and keep the level-5 table for hydrogen, carbon and nitrogen. On the
probe this reduces the bromine violation from 0.73 to below 10⁻⁶
Eh/Bohr², i.e. to the level of the light atoms, so that the raw and the
repaired spectra coincide. The repair is nevertheless kept switched on and
both spectra are printed: if a run's per-atom violation is not near 10⁻⁶,
the grid did not take effect, and the printout shows it.

Pruning (removing angular points near the nucleus) and the angular order
were shown not to matter, so the dense grid costs only the extra radial
shells: about 60 % more grid points for a five-atom complex and roughly a
third more Hessian time.

### B.7 Counting the external modes: linear complexes and the rotor test

For a linear complex only two rotations exist, so 3N−5 vibrations remain
and 5 external modes must be discarded; for a bent complex 6. The
re-analysis detects linearity geometrically (all atoms within 10⁻³ Å of a
line). The re-optimization inputs use PySCF's rotor test at the optimized
geometry and drop 5 or 6 modes accordingly. This matters because the
pre-2.14 fork's `rotation_const` did not zero the near-zero moment of
inertia of a linear molecule before classifying the rotor, so a linear
complex could be classified as non-linear at random and one genuine
bending mode was then dropped as a "rotation"; version 2.14.0 fixes this
by zeroing moments below 10⁻⁹. All HCN···HX global and halogen-bonded
minima are linear and were affected in the original logs.

### B.8 Two other failure modes seen in the original runs

- **CPHF non-convergence.** The CNEO coupled-perturbed equations are solved
  with PySCF's Krylov solver, which raises an exception when it reaches its
  iteration cap (100 in the fork). Every bent HCN···HX Hessian and several
  aug-cc-pVQZ ones died this way; the traceback went to the SLURM error
  file, and the log ends at `BEGIN HESSIAN CALCULATION`. The new inputs
  allow 300 iterations, tighten `conv_tol_cpscf` to 10⁻⁹, and retry once
  with a 0.2 Eh level shift on the preconditioner, which changes the
  convergence path but not the converged solution. On the linear HCN···HBr
  geometry the CPHF converged in 28 iterations to a residual of 2×10⁻⁷.
- **Optimizer oscillation.** The custom grms criterion of 3×10⁻⁶ Eh/Bohr
  in the original inputs lay below the gradient noise of a level-3 grid
  without grid response, and net-force and net-torque components of the
  noisy gradient made several optimizations oscillate for hundreds to
  1500 steps. The new inputs use grid response, a level-5 grid, the
  GAU_TIGHT set, and geomeTRIC's `subfrctor=2`, which projects the net
  force and torque out of every gradient.

---

## Part C. Suggested wording and references

### C.1 Paragraph for the methods section (adapt as needed)

"Harmonic frequencies were obtained from the analytic CNEO-DFT Hessian.
Because the exchange–correlation contribution to the Hessian is evaluated
on atom-centred quadrature grids whose points and weights are held fixed
during differentiation, the computed Hessian violates the translational
sum rule Σ_B ∂²E/∂R_A∂R_B = 0 by an isotropic error confined to the
diagonal block of each atom, whose magnitude is a quadrature error that
grows steeply with nuclear charge (≈10⁻³ Eh/Bohr² for F, 0.02–0.4 for Cl
and 1–18 for Br on PySCF's level-3 grid) and vanishes as the radial grid
of that atom is refined (0.73, 6.7×10⁻⁴ and 9×10⁻⁷ Eh/Bohr² for Br with
120, 240 and 400 radial shells). Left uncorrected, this error displaces the
translational modes to finite frequency and, after projection of the
external modes, contaminates the low-frequency intermolecular vibrations,
producing spurious imaginary frequencies for the Cl and Br complexes. We
therefore enforced the translational sum rule by resetting each diagonal
block to minus the sum of the corresponding off-diagonal blocks, the
standard acoustic-sum-rule correction of lattice dynamics, and verified
that the corrected spectrum reproduces that of a radially converged grid to
0.1 cm⁻¹. In the final calculations the halogen atoms carried 400 radial
shells, which removes the defect at source; the correction was retained as
a check and changed no frequency by more than [x] cm⁻¹."

### C.2 References

Verify volume and page numbers against the originals before submission.

Sum rules, Eckart conditions and projected Hessians

1. M. Born and K. Huang, *Dynamical Theory of Crystal Lattices* (Oxford University Press, 1954). Translational invariance of force constants; the acoustic sum rule.
2. E. B. Wilson, J. C. Decius and P. C. Cross, *Molecular Vibrations* (McGraw-Hill, 1955). Eckart conditions, separation of translations and rotations, normal-mode analysis.
3. W. H. Miller, N. C. Handy and J. E. Adams, "Reaction path Hamiltonian for polyatomic molecules," *J. Chem. Phys.* **72**, 99 (1980). The projected Hessian (projection of translations and rotations, and of the gradient direction).
4. X. Gonze and C. Lee, "Dynamical matrices, Born effective charges, dielectric permittivity tensors, and interatomic force constants from density-functional perturbation theory," *Phys. Rev. B* **55**, 10355 (1997). Violation of the acoustic sum rule by discretization and its enforcement.
5. N. Mounet and N. Marzari, "First-principles determination of the structural, vibrational and thermodynamic properties of diamond, graphite, and derivatives," *Phys. Rev. B* **71**, 205214 (2005). Discussion of acoustic-sum-rule imposition on computed force constants.
6. P. Giannozzi *et al.*, "QUANTUM ESPRESSO: a modular and open-source software project for quantum simulations of materials," *J. Phys.: Condens. Matter* **21**, 395502 (2009); and the `matdyn.x`/`dynmat.x` documentation, option `asr` (`'simple'` = diagonal correction, `'crystal'` = optimized projection).

Grid-weight derivatives and their effect on DFT derivatives and frequencies

7. B. G. Johnson, P. M. W. Gill and J. A. Pople, "The performance of a family of density functional methods," *J. Chem. Phys.* **98**, 5612 (1993). Grid-weight derivatives in DFT gradients.
8. B. G. Johnson and M. J. Frisch, "An implementation of analytic second derivatives of the gradient-corrected density functional energy," *J. Chem. Phys.* **100**, 7429 (1994). Analytic DFT Hessians and the role of grid terms.
9. J. Baker, J. Andzelm, A. Scheiner and B. Delley, "The effect of grid quality and weight derivatives in density functional calculations," *J. Chem. Phys.* **101**, 8894 (1994). Loss of translational/rotational invariance when weight derivatives are omitted; grid-quality dependence.
10. M. Malagoli and J. Baker, "The effect of grid quality and weight derivatives in density functional calculations of harmonic vibrational frequencies," *J. Chem. Phys.* **119**, 12763 (2003). Directly on point: frequency errors and spurious low modes from omitted weight derivatives, and their grid dependence.

Integration grids

11. A. D. Becke, "A multicenter numerical integration scheme for polyatomic molecules," *J. Chem. Phys.* **88**, 2547 (1988).
12. O. Treutler and R. Ahlrichs, "Efficient molecular numerical integration schemes," *J. Chem. Phys.* **102**, 346 (1995). Radial grids and atomic-size adjustments used by PySCF.
13. V. I. Lebedev and D. N. Laikov, "A quadrature formula for the sphere of the 131st algebraic order of accuracy," *Dokl. Math.* **59**, 477 (1999). Angular grids.

Methods and software

14. X. Xu and Y. Yang, "Constrained nuclear-electronic orbital density functional theory: Energy surfaces with nuclear quantum effects," *J. Chem. Phys.* **152**, 084107 (2020). CNEO-DFT.
15. X. Xu and Y. Yang, "Full-quantum descriptions of molecular systems from constrained nuclear–electronic orbital density functional theory," *J. Chem. Phys.* **153**, 074106 (2020). CNEO-DFT geometries and vibrational analysis.
16. X. Xu, Z. Chen and Y. Yang, "Molecular dynamics with constrained nuclear electronic orbital density functional theory: Accurate vibrational spectra from efficient incorporation of nuclear quantum effects," *J. Am. Chem. Soc.* **144**, 4039 (2022). CNEO vibrational frequencies (check title/citation).
17. Q. Yu, F. Pavošević and S. Hammes-Schiffer, "Development of nuclear basis sets for multicomponent quantum chemistry methods," *J. Chem. Phys.* **152**, 244123 (2020). The PB4-D protonic basis.
18. Q. Sun *et al.*, "Recent developments in the PySCF program package," *J. Chem. Phys.* **153**, 024109 (2020). PySCF.
19. L.-P. Wang and C. Song, "Geometry optimization made simple with translation and rotation coordinates," *J. Chem. Phys.* **144**, 214108 (2016). geomeTRIC and TRIC coordinates.
20. R. A. Kendall, T. H. Dunning Jr. and R. J. Harrison, *J. Chem. Phys.* **96**, 6796 (1992) (aug-cc-pVXZ, H and first row); D. E. Woon and T. H. Dunning Jr., *J. Chem. Phys.* **98**, 1358 (1993) (Al–Ar); A. K. Wilson, D. E. Woon, K. A. Peterson and T. H. Dunning Jr., *J. Chem. Phys.* **110**, 7667 (1999) (Ga–Kr, all-electron).
21. Functionals: J. P. Perdew, K. Burke and M. Ernzerhof, *Phys. Rev. Lett.* **77**, 3865 (1996) (PBE); J. P. Perdew *et al.*, *Phys. Rev. B* **46**, 6671 (1992) (PW91); A. D. Becke, *Phys. Rev. A* **38**, 3098 (1988) and J. P. Perdew, *Phys. Rev. B* **33**, 8822 (1986) (BP86); C. Lee, W. Yang and R. G. Parr, *Phys. Rev. B* **37**, 785 (1988) (LYP); A. D. Becke, *J. Chem. Phys.* **107**, 8554 (1997) (B97); S. Lehtola *et al.*, *SoftwareX* **7**, 1 (2018) (libxc).
