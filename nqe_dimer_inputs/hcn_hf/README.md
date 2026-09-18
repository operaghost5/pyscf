# HCN···HF — CNEO-DFT re-optimization, Hessian and harmonic frequencies

Thirty self-contained PySCF drivers: five functionals (PBE, PW91, BP86,
BLYP, B97) × three electronic bases (aug-cc-pVDZ, aug-cc-pVTZ,
aug-cc-pVQZ) × two isomers of the hydrogen cyanide / hydrogen fluoride
dimer, both hydrogen nuclei quantum with the PB4-D protonic basis, no
electron-proton correlation functional. Each input

1. runs a tight CNEO-DFT SCF at the starting geometry,
2. re-optimizes it with geomeTRIC and the analytic CNEO-DFT gradient,
3. runs a single point at the optimized geometry,
4. computes the analytic CNEO-DFT Hessian,
5. does the harmonic analysis, and
6. prints a summary: method, functional, quantum nuclei, EPC, electronic
   and nuclear basis, nuclear masses, the optimized geometry (symbols and
   Cartesians, Angstrom), the single-point energy, the number of imaginary
   modes, all 3N frequencies with the translations/rotations tagged, and
   the vibrational frequencies with those modes dropped.

```bash
python HCNHF.cneo-dft-pbe.global-min.aug-cc-pvdz.pb4d.2026-09-18.py
```

Next to the log it writes `<stem>.traj.xyz` (every optimizer step),
`<stem>.opt.xyz` and `<stem>.hessian.npz` (Hessian, frequencies,
geometry, energies).

## File naming

```
HCNHF.cneo-dft-<xc>.<global-min|local-min>.<electronic basis>.pb4d.<date>.py
```

`global-min` is the linear H–C≡N···H–F structure (HF the hydrogen-bond
donor, N the acceptor). `local-min` is the bent N≡C–H···F–H structure
(the C–H the donor, F the acceptor), which the February runs found 4.8
to 6.0 kcal/mol above the linear isomer. Atom order is H, C, N, H, F in
every file; atoms 0 and 3 are the quantum protons.

## Starting geometries

`geometries/hcn_hf_start_geometries.xyz` holds the 30 starting points:
the final geometries of the February 2026 CNEO-DFT runs of the same
functional and basis, extracted from those logs (Angstrom). The comment
line of each block records the run's final energy and whether its
optimizer converged. Three of them are not converged stationary points
and say so: B97/aug-cc-pVDZ (both isomers) and B97/aug-cc-pVTZ
(global-min) hit the 1500-step cap, and B97/aug-cc-pVQZ (global-min) is
the geometry at the last completed step of a run killed at step 358.
All are close enough to the stationary point that the re-optimization
here is a short relaxation.

## Why these convergence and grid settings

The February runs used `convergence_grms = 3e-6` Eh/Bohr with otherwise
Gaussian-default criteria, the default level-3 grid, the default SCF
tolerance (`conv_tol = 1e-9`, which sets the orbital-gradient threshold
to 3e-5), and a 1500-step cap. Decomposing the final gradients of those
runs showed that the residual "gradient" was dominated by a net force
and torque of 7e-6 to 5e-5 Eh/Bohr, which is a numerical artefact (grid
translational-invariance error plus incomplete SCF) that no geometry
change can remove. The linear B97 runs chased it for 1500 steps; the
bent runs only converged because geomeTRIC's auto-detection started
projecting the net force and torque out. The settings in the inputs fix
both the noise and the test:

| Setting | February 2026 | Here | Purpose |
|---|---|---|---|
| geomeTRIC criteria | grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3, ΔE 1e-9 | `GAU_TIGHT` (grms 1e-5, gmax 1.5e-5, drms 4e-5, dmax 6e-5) with ΔE 1e-8 | self-consistent set that is adequate for the 80 cm⁻¹ modes |
| `subfrctor` | 1 (auto) | 2 (always) | project net force and torque out of the gradient |
| optimizer step cap | 1500 | 200 | a 5-atom complex not done in 200 steps is stuck in noise |
| SCF `conv_tol` / `conv_tol_grad` | 1e-9 / 3.2e-5 | 1e-11 / 1e-6 | gradient error is first order in the orbital-gradient residual |
| SCF `max_cycle` | 50 | 200 | headroom for the tighter SCF (an unconverged SCF aborts the optimizer) |
| electronic grid level | 3 | 5 (`mf.components['e'].grids.level`) | smaller grid error in energy, gradient and Hessian |
| `grid_response` | off | on | removes the grid's net-force error; allowed because EPC is off |
| `conv_tol_cpscf` | 1e-8 | 1e-9 | tighter CPHF for the Hessian's low-frequency modes |
| Hessian CPHF `max_cycle` | 100 | 300, then one retry with `level_shift` 0.2 Eh | see below |

**Why the February Hessians failed.** Seventeen of the February runs
stopped inside the Hessian. Every one of them stopped after exactly 300
inter-component potential builds, three per iteration of the
coupled-perturbed solver, i.e. at the solver's 100-iteration cap, while
the ten Hessians that finished needed 59 to 69 iterations. PySCF's
Krylov solver raises `RuntimeError("Krylov solver failed to converge.")`
at the cap, so those runs crashed, with the traceback in the `.err`
file, rather than running out of walltime. The failures were every
aug-cc-pVQZ Hessian and every bent-isomer Hessian, so CNEO-CPHF
convergence degrades with basis size and for the bent structure. The
inputs therefore allow 300 iterations, log the Krylov residual each
iteration (`run_hessian.verbose = 4`), and on failure retry once with a
0.2 Eh level shift on the preconditioner, which does not change the
converged result.

Every setting is a named constant at the top of each input, and
`STOP_IF_OPT_FAILS = True` stops a run before the Hessian if geomeTRIC
gives up, instead of spending the Hessian on an unconverged geometry.

## Frequency analysis

The Hessian is analysed twice with the fork's `Hessian.harmonic_analysis`,
which uses the CNEO nuclear masses (`mol.mass`):

- **all 3N modes**, nothing projected, ascending, negative = imaginary,
  with the smallest |f| modes tagged `TR`; then the same list with the
  TR modes dropped. The number of TR modes is **5 for a linear rotor and
  6 otherwise**, decided with PySCF's own rotor test at the optimized
  geometry. Dropping six for the linear isomer would discard one of its
  two degenerate 80 cm⁻¹ bends.
- **PySCF's projected analysis** (translations and rotations removed
  before diagonalization) as a one-line cross-check; the two lists
  should agree to within the size of the TR modes.

The translational sum rule of the Hessian is printed as a quality check;
the external modes of the February Hessians were 20 to 40 cm⁻¹, and the
tighter SCF/grid here should shrink both.

## Notes

- PySCF maps `xc='b97'` to libxc's `HYB_GGA_XC_B97`, Becke's 1997 hybrid
  with about 19 % exact exchange, not a pure GGA. This matches the
  February runs and explains B97's systematically higher stretching
  frequencies.
- The nuclear masses printed in the summary are the CNEO values: most
  common isotope minus the electron mass for the quantum protons,
  isotope-averaged for the classical nuclei. The February frequency
  analyses used isotope-averaged masses for every atom.
- Regenerate the inputs with `python make_inputs.py [date]` after editing
  the template in `make_inputs.py`; do not edit the inputs by hand.

## Running on the cluster

```bash
./submit_array.slurm
```

on the login node submits three arrays (dz/tz/qz tiers, resources in the
script header, all within a 2-day queue limit) over the files in
`inputs/`; `INPUT_DIR`, `OUTPUT_DIR` and `PYTHON` are overridable in the
environment. Each task checks that `pyscf.neo` and `geometric` import and
warns if the PySCF build lacks the linear-rotor fix.

Cost guide from the February logs, which ran single-threaded: at
aug-cc-pVQZ one SCF took about 200 s and one gradient about 340 s, so
the re-optimization plus single point is 1 to 3 hours there even with
the finer grid, and the Hessian is dominated by the coupled-perturbed
solve at roughly 15 Fock-like builds per Krylov iteration. On the
32-thread QZ tier the whole input should fit comfortably inside 47
hours provided the CPHF converges; the DZ inputs take minutes and the TZ
inputs about an hour.
