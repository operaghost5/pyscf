#!/usr/bin/env python
"""
Probe: why did every HCN...HBr aug-cc-pVTZ Hessian die at 'BEGIN HESSIAN CALCULATION'?

Runs CNEO-DFT/PBE aug-cc-pVTZ at the February 2026 optimized linear geometry,
then executes the Hessian STAGE BY STAGE in the same order as
pyscf.neo.hessian.hess_elec, printing a marker before and after each stage so
the failing stage is identified even without the SLURM .err file:

  A  Hessian object construction
  B  partial_hess_int   (inter-component 2nd-derivative integrals, the first
                         thing hess_elec does; the February logs died before it
                         printed anything)
  C  make_h1            (first-derivative Fock matrices)
  D  solve_mo1          (CNEO-CPHF, Krylov)
  E  per-component hess_elec (electronic ERI + XC second derivatives)
  F  hess_nuc, assembly, harmonic analysis

faulthandler is enabled on stdout, so a hard crash (segmentation fault, illegal
instruction, abort) also leaves a Python traceback in the .out. Python
exceptions are caught and printed with their traceback.

Run with the same environment as the production jobs, unbuffered:

    ulimit -s unlimited          # PySCF/libcint recommendation; try WITHOUT it first
    python -u hcnhbr_avtz_hessian_probe.py > probe.out 2>&1

Expected wall time on one node: SCF ~30 s, stage B a few minutes, stage D up to
an hour if the CPHF runs to 300 iterations.
"""
import faulthandler
import os
import platform
import resource
import sys
import time
import traceback

faulthandler.enable(file=sys.stdout, all_threads=True)

import numpy
import pyscf
from pyscf import neo

XC = 'pbe'
BASIS = 'aug-cc-pvtz'
NUC_BASIS = 'pb4d'
MAX_MEMORY = int(os.environ.get('PYSCF_MAX_MEMORY', 16000))
HESS_MAX_CYCLE = 300

# February 2026 PBE/aug-cc-pVTZ optimized linear H-C#N...H-Br geometry (Angstrom)
GEOM = """
H   0.0  0.0   4.66585492
C   0.0  0.0   3.56547366
N   0.0  0.0   2.40962255
H   0.0  0.0   0.46672139
Br  0.0  0.0  -1.01691374
"""


def stamp(msg):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    print(f'[{time.strftime("%H:%M:%S")}] {msg}   (max RSS so far {rss:.0f} MB)', flush=True)


def cpu_info():
    model = flags = 'n/a'
    try:
        for line in open('/proc/cpuinfo'):
            if line.startswith('model name') and model == 'n/a':
                model = line.split(':', 1)[1].strip()
            if line.startswith('flags') and flags == 'n/a':
                fl = line.split(':', 1)[1].split()
                flags = ' '.join(f for f in fl if f.startswith(('avx', 'sse4', 'fma')))
    except OSError:
        pass
    return model, flags


def run_stage(label, func):
    stamp(f'stage {label}: start')
    try:
        out = func()
    except BaseException:
        stamp(f'stage {label}: FAILED with a Python exception')
        traceback.print_exc(file=sys.stdout)
        sys.stdout.flush()
        sys.exit(1)
    stamp(f'stage {label}: done')
    return out


# ----------------------------------------------------------------------------
print('probe: HCN...HBr CNEO-DFT/PBE aug-cc-pVTZ Hessian, stage by stage', flush=True)
print(f'host {platform.node()}   python {sys.version.split()[0]}   numpy {numpy.__version__}')
print(f'pyscf {pyscf.__version__} from {os.path.dirname(pyscf.__file__)}')
model, flags = cpu_info()
print(f'cpu {model}\n    simd flags: {flags}')
soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
print(f'stack limit soft/hard: {soft if soft != resource.RLIM_INFINITY else "unlimited"} / '
      f'{hard if hard != resource.RLIM_INFINITY else "unlimited"}')
print(f'OMP_NUM_THREADS={os.environ.get("OMP_NUM_THREADS")}   pyscf threads {pyscf.lib.num_threads()}   '
      f'max_memory {MAX_MEMORY} MB', flush=True)

mol = neo.M(atom=GEOM, basis=BASIS, nuc_basis=NUC_BASIS, quantum_nuc=['H'],
            charge=0, unit='Angstrom', verbose=4, max_memory=MAX_MEMORY)
print(f'electronic AOs: {mol.components["e"].nao}   electrons: {mol.components["e"].nelectron}   '
      f'quantum nuclei: {mol.nuc_num}', flush=True)

mf = neo.CDFT(mol, xc=XC, epc=None)
mf.conv_tol = 1e-11
mf.conv_tol_grad = 1e-6
mf.max_cycle = 200
mf.conv_tol_cpscf = 1e-9
mf.components['e'].grids.level = 5
stamp('SCF: start')
e_scf = mf.kernel()
stamp(f'SCF: done, E = {e_scf:.10f} Eh, converged = {mf.converged}')

mo_energy, mo_coeff, mo_occ = mf.mo_energy, mf.mo_coeff, mf.mo_occ
atmlst = list(range(mol.natm))


def stage_a():
    h = mf.Hessian()
    h.max_cycle = HESS_MAX_CYCLE
    h.verbose = 5
    return h


hess = run_stage('A (Hessian object)', stage_a)
log = pyscf.lib.logger.new_logger(hess, 5)

de2 = run_stage('B (partial_hess_int: inter-component derivative integrals)',
                lambda: hess.partial_hess_int(mo_coeff, mo_occ, atmlst, log))
h1ao = run_stage('C (make_h1: first-derivative Fock matrices)',
                 lambda: hess.make_h1(mo_coeff, mo_occ, None, atmlst, log))
mo1, mo_e1 = run_stage('D (solve_mo1: CNEO-CPHF)',
                       lambda: hess.solve_mo1(mo_energy, mo_coeff, mo_occ, h1ao,
                                              None, atmlst, MAX_MEMORY, log))
for t, comp in hess.components.items():
    de2 = de2 + run_stage(f'E (hess_elec of component {t})',
                          lambda t=t, comp=comp: comp.hess_elec(
                              mo_energy[t], mo_coeff[t], mo_occ[t], mo1[t], mo_e1[t],
                              h1ao[t], atmlst, MAX_MEMORY, log))
de = run_stage('F (hess_nuc + assembly)', lambda: de2 + hess.hess_nuc(mol, atmlst))

viol = numpy.abs(de.sum(axis=1)).reshape(mol.natm, -1).max(axis=1)
print('translational sum rule per atom (Eh/Bohr^2): '
      + ', '.join(f'{mol.atom_pure_symbol(i)}{i} {v:.2e}' for i, v in enumerate(viol)))
res = hess.harmonic_analysis(mol, de, exclude_trans=False, exclude_rot=False,
                             imaginary_freq=False, intensity=False)
freqs = numpy.sort(numpy.asarray(res['freq_wavenumber'], dtype=float))
print('all 3N frequencies, raw Hessian (cm^-1): '
      + ', '.join(f'{f:.1f}' for f in freqs))
stamp('ALL STAGES OK')
