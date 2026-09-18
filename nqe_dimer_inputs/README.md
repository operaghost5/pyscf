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

## Requirements

- `theorychemyang/pyscf` at version 2.14 or later, built from source.
  Earlier builds mis-classify a linear molecule as non-linear at random
  in `pyscf.hessian.thermo` and drop one real bending mode.
- geomeTRIC (any 1.x).
- Do not patch library files for run settings; every setting the inputs
  need is exposed on the PySCF objects and set in the driver scripts.
