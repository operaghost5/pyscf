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

## Requirements

- `theorychemyang/pyscf` at version 2.14 or later, built from source.
  Earlier builds mis-classify a linear molecule as non-linear at random
  in `pyscf.hessian.thermo` and drop one real bending mode.
- geomeTRIC (any 1.x).
- Do not patch library files for run settings; every setting the inputs
  need is exposed on the PySCF objects and set in the driver scripts.
