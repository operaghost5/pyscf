# HCN···HBr aug-cc-pVTZ Hessian probe — results (2026-09-18)

Question: why did every February 2026 HCN···HBr aug-cc-pVTZ Hessian die at
`BEGIN HESSIAN CALCULATION` with nothing after it and no SLURM error file?

Setup: CNEO-DFT/PBE, aug-cc-pVTZ (all-electron Br, 197 AOs, 50 electrons),
PB4-D on both protons, no EPC, at the February PBE/aug-cc-pVTZ linear
geometry; SCF `conv_tol` 1e-11, `conv_tol_grad` 1e-6, `conv_tol_cpscf` 1e-9,
grid level 5 (176 912 points), Hessian `max_cycle` 300. PySCF fork 2.14.0,
Python 3.12.14, numpy 2.5.2, Xeon Gold 6138 nodes.

## Run 1: one core, interactive (`hcnhbr_avtz_hessian_probe.py`)

| Stage | Wall | Peak RSS |
|---|---|---|
| SCF, 14 cycles | 37 s | 2.2 GB |
| B inter-component second-derivative integrals | 19 s | 2.2 GB |
| C first-derivative Fock matrices | 102 s | 3.7 GB |
| D CNEO-CPHF, 28 Krylov cycles, residual 2e-7 | 530 s | 3.7 GB |
| E electronic second derivatives | 1123 s | 4.0 GB |
| whole Hessian | 30 min | 4.0 GB |

Sum rule, Eh/Bohr²: H 7e-8, C 9e-6, N 1.2e-5, H 1.2e-7, Br 7.3e-1.
Raw 3N spectrum, cm⁻¹: -1.8, -1.8, 43.2, 43.2, 91.2, 437.8, 437.8, 491.4,
494.8, 494.8, 710.3, 710.3, 2129.4, 2223.8, 3217.2. The three translations
sit at 438/491 cm⁻¹; an isotropic error δ on the Br block moves them to
about sqrt(δ/M_total) = 423 cm⁻¹, which confirms the diagonal-block model.

Conclusion: neither the Hessian's time nor its memory can have killed the
February jobs on any reasonable allocation. The footprint steps from 2.2 GB
(optimization) to 3.7 GB (Hessian start), so a memory limit between those
two numbers produces exactly the observed silent death; walltime spent by
a long optimization and preemption on the `requeue` partition are the
other candidates. `sacct` State/ReqMem/Elapsed decide between them.

## Run 2: 16 cores, batch, job 17812751 (`hcnhbr_avtz_hessian_probe2.py`)

cgroup limit 68 719 MB, PYSCF_MAX_MEMORY 32 000.

| Stage | Wall | Peak RSS |
|---|---|---|
| SCF | 7 s | 2.2 GB |
| B | 8 s | 2.2 GB |
| C | 17 s | 3.8 GB |
| D, 28 Krylov cycles | 2.1 min | 3.8 GB |
| E, of which ERIs 2.0 min, XC diag 2 s, XC deriv2 15 s | 2.3 min | 4.1 GB |
| whole Hessian | 4.9 min | 4.1 GB |

Speed-up over one core: 6× overall (ERIs 8.7×, CPHF 4.5×).

### Stage G: bromine grid test (SCF + full Hessian each)

| Grid | points | Br violation | raw 3N spectrum, cm⁻¹ |
|---|---|---|---|
| level 5 (Br 120 × 770) | 176 912 | 7.32e-1 | -1.8, -1.8, 43.2, 43.2, 91.2, 437.8, 437.8, 491.4, 494.8, 494.8, 710.3, 710.3, 2129.4, 2223.8, 3217.2 |
| level 5, no pruning | 336 704 | 7.32e-1 | identical |
| Br 240 × 770 | 221 104 | 6.68e-4 | -1.8, -1.8, 12.8, 14.7, 14.7, 43.4, 43.4, 106.2, 441.0, 441.0, 710.3, 710.3, 2129.4, 2223.0, 3217.2 |
| Br 240 × 770, no pruning | 429 104 | 6.68e-4 | identical |
| Br 400 × 974, no pruning | 633 904 | 8.98e-7 | -1.9, -1.9, 0.7, 0.7, 0.8, 43.4, 43.4, 105.9, 441.0, 441.0, 710.3, 710.3, 2129.4, 2223.0, 3217.2 |

The other atoms' violations are unchanged (C 9e-6, N 1.2e-5, H ≤ 1e-7),
the SCF energies agree to 0.03 µEh, and the repaired spectrum is the same
for every grid: 0.0, 0.0, 0.0, 0.4, 0.4, 43.3, 43.3, 105.9, 441.0, 441.0,
710.3, 710.3, 2129.4, 2223.0, 3217.2 cm⁻¹.

Conclusions:

1. The bromine defect is radial quadrature stiffness of the Br core on a
   grid that the Hessian formulas hold fixed. Pruning and the angular
   order are irrelevant; 120 → 240 → 400 radial shells take it from 0.73
   to 6.7e-4 to 9e-7 Eh/Bohr².
2. The diagonal-block repair is exact for this defect: the repaired
   120-shell spectrum equals the raw 400-shell spectrum to 0.1 cm⁻¹ in
   every vibration.
3. Assignment of the ten vibrations (linear complex, 5 TR modes):
   intermolecular bend 43 (π), intermolecular stretch 106 (σ), H–Br
   libration 441 (π), HCN bend 710 (π), C≡N stretch 2129, H–Br stretch
   2223, C–H stretch 3217 cm⁻¹. On the uncorrected level-5 Hessian the
   intermolecular stretch appeared at 91 and the libration at 495 cm⁻¹
   through mixing with the displaced translations.
4. The production inputs use `BR_ATOM_GRID = (400, 770)` with default
   pruning and keep the repair as a check.
