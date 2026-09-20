# NQE dimer inputs

Nuclear quantum effects in hydrogen-bonded dimers with constrained
nuclear-electronic orbital DFT (CNEO-DFT, `pyscf.neo.CDFT` from the
`theorychemyang/pyscf` fork). One folder per system; each holds the
starting geometries, a generator that writes the PySCF driver scripts,
the generated inputs, a SLURM array script and a README with the
protocol and its provenance.

## Contents

- `hcn_hf/` — HCN···HF. Re-optimization, analytic Hessian and harmonic
  frequencies of the linear N···H–F global minimum and the bent
  C–H···F–H local minimum with PBE, PW91, BP86, BLYP and B97 in
  aug-cc-pVDZ/TZ/QZ, both protons quantum (PB4-D), starting from the
  February 2026 CNEO-DFT geometries. See `hcn_hf/README.md`.
- `hcn_hcl/` — HCN···HCl. The same protocol for the linear N···H–Cl
  global minimum, the bent C–H···Cl–H local minimum and the linear
  halogen-bonded N···Cl–H isomer, plus a translational sum-rule repair of
  the Hessian that removes a chlorine-localized defect responsible for
  spurious imaginary bends in the February runs. See `hcn_hcl/README.md`.
- `hcn_hbr/` — HCN···HBr. The same three isomers with bromine
  all-electron, plus a 400-shell radial grid on bromine: a stage-by-stage
  probe (`hcn_hbr/probe/`) showed the bromine sum-rule defect of the
  fixed-grid Hessian is radial quadrature error (0.73 → 6.7e-4 → 9e-7
  Eh/Bohr² for 120 → 240 → 400 shells) and that the diagonal-block repair
  reproduces the converged spectrum to 0.1 cm⁻¹. See `hcn_hbr/README.md`.
- `hcn_hcl_densegrid/` — HCN···HCl again with 400 radial shells on
  chlorine (`cl400` file-name token, own input/output directories and
  job names); otherwise identical to `hcn_hcl/`, which is unchanged.
  See `hcn_hcl_densegrid/README.md`.
- `hf_hf/` — (HF)₂. The bent C_s minimum re-optimized from the 15 April
  2026 geometries plus a constructed C₂h interchange-saddle start, both
  protons quantum (April had only the donor), 400 radial shells on
  fluorine, grid details in the printed summary. `hf_hf/april_analysis/`
  holds the classification of the April runs (15 minima; the 20
  local-minimum runs died on an undefined geometry). See `hf_hf/README.md`.
- `hf_hbr/` — HF/HBr. The two isomers F–H···Br–H (`hf-donor`) and
  Br–H···F–H (`hbr-donor`) re-optimized from the 30 April 2026 geometries
  with both protons quantum (April had only each isomer's donor proton, so
  its energies cannot rank them), 400 radial shells on Br and F, grid
  details in the summary. `hf_hbr/april_analysis/`: all 30 April runs are
  minima once the bromine defect is repaired; the raw aug-cc-pVQZ Hessians
  showed a spurious ~400i mode. See `hf_hbr/README.md`.
- `hf_hcl/` — HF/HCl. The two isomers F–H···Cl–H (`hf-donor`) and
  Cl–H···F–H (`hcl-donor`) re-optimized from the 30 April 2026 geometries
  with both protons quantum (April had only the HCl proton quantum in
  both isomers), one atom order for both isomers, 400 radial shells on Cl
  and F, grid details in the summary. `hf_hcl/april_analysis/`: all 30
  April runs are minima once the chlorine defect is repaired; raw, the 20
  aug-cc-pVDZ and aug-cc-pVQZ runs showed a spurious 54i to 307i mode.
  See `hf_hcl/README.md`.
- `hcl_hbr/` — HCl/HBr. The two isomers Cl–H···Br–H (`hcl-donor`) and
  Br–H···Cl–H (`hbr-donor`) re-optimized from the April 2026 geometries
  with both protons quantum (April had only each isomer's donor proton, so
  its energies cannot rank them), 400 radial shells on Cl and Br, grid
  details in the summary. `hcl_hbr/april_analysis/`: all 29 finished April
  runs are minima once the bromine defect is repaired; raw, the nine
  aug-cc-pVQZ Hessians showed a spurious ~700i mode. See `hcl_hbr/README.md`.
- `hcl_hcl/` — (HCl)₂. The bent C_s minimum re-optimized from the 15 April
  2026 geometries plus a constructed C₂h interchange-saddle start, both
  protons quantum (April had only the acceptor proton), 400 radial shells
  on chlorine, grid details in the summary. `hcl_hcl/april_analysis/`: all
  15 April runs are minima once the chlorine defect is repaired (raw, the
  DZ and QZ runs showed a spurious 206i to 555i mode); the 20 local-minimum
  runs died on an undefined geometry. See `hcl_hcl/README.md`.
- `hbr_hbr/` — (HBr)₂. The bent C_s minimum re-optimized from the 15 April
  2026 geometries plus a constructed C₂h interchange-saddle start, both
  protons quantum (April had only the acceptor proton), 400 radial shells
  on bromine, grid details in the summary. `hbr_hbr/april_analysis/`: all
  15 April runs are minima once the bromine defect is repaired; raw, every
  spectrum was wrong (spurious 690 to 2327 cm⁻¹ modes at DZ/TZ, 934i to
  1045i at QZ); the 20 local-minimum runs died on an undefined geometry.
  See `hbr_hbr/README.md`.

- `hcn_hf_dft/` — HCN···HF with conventional DFT (all nuclei classical),
  the reference for the nuclear quantum effect. The same 30 linear and
  bent starts re-optimized from the February 2026 DFT geometries with
  every setting of `hcn_hf/` applied to `dft.RKS`, the CNEO mass
  convention in the harmonic analysis, the sum-rule check and repair,
  and the grid printed in the summary. `hcn_hf_dft/february_analysis/`:
  all 30 February DFT runs are minima (level-9 grid, fluorine defect
  below 3e-4 Eh/Bohr², yet the softest bent-complex mode moved by up to
  26 cm⁻¹ under the repair), plus CNEO-DFT minus DFT shifts for the ten
  runs with Hessians in both sets. See `hcn_hf_dft/README.md`.
- `METHODS.md` — computational details for every system and both
  generations of runs (software, model, quantum nuclei, grids, optimizer
  and SCF/CPHF criteria, analysis protocol, sum-rule violations found),
  a graduate-level explanation of the translational sum-rule defect of
  the fixed-grid Hessian and its diagonal-block repair, suggested methods
  wording and references.

## Requirements

- `theorychemyang/pyscf` at version 2.14 or later, built from source.
  Earlier builds mis-classify a linear molecule as non-linear at random
  in `pyscf.hessian.thermo` and drop one real bending mode.
- geomeTRIC (any 1.x).
- Do not patch library files for run settings; every setting the inputs
  need is exposed on the PySCF objects and set in the driver scripts.
