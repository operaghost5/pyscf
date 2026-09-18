#!/usr/bin/env python
"""
Probe 2: why did every HCN...HBr aug-cc-pVTZ Hessian die at 'BEGIN HESSIAN CALCULATION'?

Same calculation as the first probe (CNEO-DFT/PBE aug-cc-pVTZ at the February
2026 linear geometry, Hessian executed stage by stage), with four additions so
that the run leaves a usable record however it ends:

  * a heartbeat thread prints elapsed time, current and peak RSS, and the
    cgroup memory usage against the job's cgroup limit every few minutes.  A
    kill by the cgroup OOM killer is a SIGKILL and can leave no traceback, but
    the last heartbeat shows how close to the limit the job was and in which
    routine;
  * SIGTERM (walltime kill) and SIGUSR1 (SLURM's --signal warning) dump the
    traceback of every thread plus the memory state, so a walltime kill shows
    exactly which integral or grid routine was running;
  * stage E of the electronic component is split by wrapping the PySCF
    routines it calls, in execution order:
        E.ejk        second-derivative ERIs (int2e_ipip1, then int2e_ip1ip2
                     and int2e_ipvip1 per atom), direct contraction
        E.vxc_diag   XC second derivatives, diagonal (same-atom) blocks
        E.vxc_deriv2 XC second derivatives on the grid, all atom pairs
        (rest of E) contraction with the CPHF solution mo1
    every direct ERI contraction (get_jk) is stamped as well;
  * the Hessian logger runs at verbose 6 so PySCF's own per-atom integral
    timers and the Krylov residual lines print;
  * a table of wall time and peak RSS per stage at the end (or at the point of
    a Python exception).

Stages, as before:
  A  Hessian object construction
  B  partial_hess_int   (inter-component 2nd-derivative integrals)
  C  make_h1            (first-derivative Fock matrices)
  D  solve_mo1          (CNEO-CPHF, Krylov)
  E  per-component hess_elec (electronic ERI + XC second derivatives, then the
                              two quantum protons)
  F  hess_nuc, assembly, sum rule, harmonic analysis (raw and repaired)
  G  grid test for the bromine sum-rule defect: SCF + full Hessian again with
     pruning off, with a denser Br radial grid, and with both (see GRID_VARIANTS)

The first probe completed on one core in 31 minutes with a peak RSS of 4.0 GB
(stage D 9 min, stage E 19 min), so this run is about the multi-core timings,
the repaired spectrum and stage G.  Run it as a batch job with run_probe.sbatch
(47 h, 64 GB, 16 cores).  Environment knobs:

    PYSCF_MAX_MEMORY      MB handed to PySCF (default 16000; sbatch sets 32000)
    PROBE_HEARTBEAT_MIN   minutes between heartbeats (default 5)
"""
import faulthandler
import os
import platform
import resource
import signal
import sys
import threading
import time
import traceback

faulthandler.enable(file=sys.stdout, all_threads=True)

T0 = time.time()
HEARTBEAT_MINUTES = float(os.environ.get('PROBE_HEARTBEAT_MIN', 5))
_stage = ['startup']          # name of the routine currently running
_summary = []                 # (label, minutes, peak RSS MB)


# --- memory bookkeeping ------------------------------------------------------
def rss_mb():
    """Current and peak resident set size in MB (from /proc, else getrusage)."""
    cur = peak = float('nan')
    try:
        for line in open('/proc/self/status'):
            if line.startswith('VmRSS:'):
                cur = int(line.split()[1]) / 1024.0
            elif line.startswith('VmHWM:'):
                peak = int(line.split()[1]) / 1024.0
    except OSError:
        pass
    if peak != peak:  # NaN
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    return cur, peak


def _cgroup_files():
    """(usage file, limit file) of the memory cgroup this process runs in.

    Walks up the hierarchy until a real limit is found, because SLURM sets the
    limit on the job cgroup while the step/task cgroups below it say 'max'.
    Returns None when nothing usable is found.
    """
    try:
        lines = open('/proc/self/cgroup').read().splitlines()
    except OSError:
        return None
    for line in lines:
        parts = line.split(':', 2)
        if len(parts) != 3:
            continue
        hier, ctrls, path = parts
        if hier == '0':                                   # cgroup v2
            root, usage, limit, nolimit = '/sys/fs/cgroup', 'memory.current', 'memory.max', ('max',)
        elif 'memory' in ctrls.split(','):                # cgroup v1
            root, usage, limit, nolimit = ('/sys/fs/cgroup/memory', 'memory.usage_in_bytes',
                                           'memory.limit_in_bytes', ('9223372036854771712',))
        else:
            continue
        base = os.path.normpath(root + path)
        if not os.path.exists(os.path.join(base, usage)):
            continue
        usage_file = os.path.join(base, usage)
        node = base
        while node.startswith(root) and len(node) >= len(root):
            lim = os.path.join(node, limit)
            try:
                val = open(lim).read().strip()
            except OSError:
                val = None
            if val is not None and val not in nolimit:
                return usage_file, lim
            if node == root:
                break
            node = os.path.dirname(node)
        return usage_file, None
    return None


_CG = _cgroup_files()


def cgroup_mb():
    """(usage MB, limit MB or inf) of the job's memory cgroup, or None."""
    if _CG is None:
        return None
    usage_file, limit_file = _CG
    try:
        usage = int(open(usage_file).read().strip()) / 1e6
        limit = float('inf')
        if limit_file is not None:
            limit = int(open(limit_file).read().strip()) / 1e6
        return usage, limit
    except (OSError, ValueError):
        return None


def mem_line():
    cur, peak = rss_mb()
    s = f'RSS {cur:.0f} MB, peak {peak:.0f} MB'
    cg = cgroup_mb()
    if cg is not None:
        usage, limit = cg
        lim = 'no limit' if limit == float('inf') else f'{limit:.0f} MB limit'
        s += f', cgroup {usage:.0f} MB ({lim})'
    return s


def stamp(msg):
    mins = (time.time() - T0) / 60.0
    print(f'[{time.strftime("%H:%M:%S")} +{mins:7.1f} min] {msg}   ({mem_line()})', flush=True)


def print_summary():
    print('\nstage summary (wall minutes, peak RSS MB at the end of the stage):')
    for label, mins, peak in _summary:
        print(f'  {mins:8.1f} min   {peak:8.0f} MB   {label}')
    sys.stdout.flush()


# --- heartbeat and signals ---------------------------------------------------
def _heartbeat():
    while True:
        time.sleep(HEARTBEAT_MINUTES * 60.0)
        stamp(f'heartbeat during {_stage[0]}')


threading.Thread(target=_heartbeat, daemon=True, name='heartbeat').start()


def _on_signal(signum, frame):
    name = signal.Signals(signum).name
    stamp(f'received {name} during {_stage[0]}')
    if signum == signal.SIGTERM:
        print_summary()
        os._exit(143)


for _sig in (signal.SIGTERM, signal.SIGUSR1):
    signal.signal(_sig, _on_signal)
    # C-level dump of every thread's traceback the moment the signal arrives,
    # then the Python handler above (chain=True) once the interpreter regains
    # control.  Without this a walltime kill during a long libcint call leaves
    # nothing.
    faulthandler.register(_sig, file=sys.stdout, all_threads=True, chain=True)


# --- PySCF ---------------------------------------------------------------------
import numpy                                   # noqa: E402
import pyscf                                   # noqa: E402
from pyscf import neo                          # noqa: E402
from pyscf.hessian import rhf as rhf_hess      # noqa: E402
from pyscf.hessian import rks as rks_hess      # noqa: E402


def _wrap(module, name, label=None):
    """Replace module.name with a stamped wrapper (same call, same return)."""
    orig = getattr(module, name)

    def wrapper(*args, **kwargs):
        lab = label
        if lab is None:   # _get_jk: name the integral being contracted
            lab = f'get_jk {args[1]}' if len(args) > 1 else f'{name}'
        prev = _stage[0]
        _stage[0] = lab
        stamp(f'    {lab}: start')
        t = time.time()
        out = orig(*args, **kwargs)
        stamp(f'    {lab}: done in {(time.time() - t) / 60.0:.1f} min')
        _stage[0] = prev
        return out

    wrapper.__wrapped__ = orig
    setattr(module, name, wrapper)


_wrap(rhf_hess, '_partial_hess_ejk', 'E.ejk (second-derivative ERIs, direct)')
_wrap(rks_hess, '_get_vxc_diag', 'E.vxc_diag (XC second derivatives, diagonal blocks)')
_wrap(rks_hess, '_get_vxc_deriv2', 'E.vxc_deriv2 (XC second derivatives on the grid)')
_wrap(rhf_hess, '_get_jk')

XC = 'pbe'
BASIS = 'aug-cc-pvtz'
NUC_BASIS = 'pb4d'
MAX_MEMORY = int(os.environ.get('PYSCF_MAX_MEMORY', 16000))
HESS_MAX_CYCLE = 300
GRID_LEVEL = 5

# February 2026 PBE/aug-cc-pVTZ optimized linear H-C#N...H-Br geometry (Angstrom)
GEOM = """
H   0.0  0.0   4.66585492
C   0.0  0.0   3.56547366
N   0.0  0.0   2.40962255
H   0.0  0.0   0.46672139
Br  0.0  0.0  -1.01691374
"""


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
    prev = _stage[0]
    _stage[0] = label
    stamp(f'stage {label}: start')
    t = time.time()
    try:
        out = func()
    except BaseException:
        stamp(f'stage {label}: FAILED with a Python exception')
        traceback.print_exc(file=sys.stdout)
        _summary.append((label + ' (FAILED)', (time.time() - t) / 60.0, rss_mb()[1]))
        print_summary()
        sys.exit(1)
    mins = (time.time() - t) / 60.0
    _summary.append((label, mins, rss_mb()[1]))
    stamp(f'stage {label}: done in {mins:.1f} min')
    _stage[0] = prev
    return out


# ----------------------------------------------------------------------------
print('probe 2: HCN...HBr CNEO-DFT/PBE aug-cc-pVTZ Hessian, stage by stage', flush=True)
print(f'host {platform.node()}   python {sys.version.split()[0]}   numpy {numpy.__version__}')
print(f'pyscf {pyscf.__version__} from {os.path.dirname(pyscf.__file__)}')
model, flags = cpu_info()
print(f'cpu {model}\n    simd flags: {flags}')
soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
print(f'stack limit soft/hard: {soft if soft != resource.RLIM_INFINITY else "unlimited"} / '
      f'{hard if hard != resource.RLIM_INFINITY else "unlimited"}')
slurm = {k: os.environ[k] for k in sorted(os.environ)
         if k in ('SLURM_JOB_ID', 'SLURM_JOB_PARTITION', 'SLURM_CPUS_PER_TASK', 'SLURM_CPUS_ON_NODE',
                  'SLURM_MEM_PER_NODE', 'SLURM_MEM_PER_CPU', 'SLURM_JOB_NODELIST')}
print('slurm: ' + (', '.join(f'{k}={v}' for k, v in slurm.items()) or 'not under SLURM'))
cg = cgroup_mb()
if cg is None:
    print('cgroup memory limit: not readable')
else:
    print(f'cgroup memory limit: '
          f'{"none" if cg[1] == float("inf") else f"{cg[1]:.0f} MB"}   (usage now {cg[0]:.0f} MB)')
print(f'OMP_NUM_THREADS={os.environ.get("OMP_NUM_THREADS")}   pyscf threads {pyscf.lib.num_threads()}   '
      f'max_memory {MAX_MEMORY} MB   heartbeat every {HEARTBEAT_MINUTES:g} min', flush=True)

mol = neo.M(atom=GEOM, basis=BASIS, nuc_basis=NUC_BASIS, quantum_nuc=['H'],
            charge=0, unit='Angstrom', verbose=4, max_memory=MAX_MEMORY)
print(f'electronic AOs: {mol.components["e"].nao}   electrons: {mol.components["e"].nelectron}   '
      f'quantum nuclei: {mol.nuc_num}   '
      f'nuclear AOs: {", ".join(str(mol.components[k].nao) for k in mol.components if k != "e")}',
      flush=True)

mf = neo.CDFT(mol, xc=XC, epc=None)
mf.conv_tol = 1e-11
mf.conv_tol_grad = 1e-6
mf.max_cycle = 200
mf.conv_tol_cpscf = 1e-9
mf.components['e'].grids.level = GRID_LEVEL
e_scf = run_stage('SCF', mf.kernel)
print(f'SCF: E = {e_scf:.10f} Eh, converged = {mf.converged}, '
      f'grid points {mf.components["e"].grids.weights.size}', flush=True)

mo_energy, mo_coeff, mo_occ = mf.mo_energy, mf.mo_coeff, mf.mo_occ
atmlst = list(range(mol.natm))


def stage_a():
    h = mf.Hessian()
    h.max_cycle = HESS_MAX_CYCLE
    h.verbose = 6
    return h


hess = run_stage('A (Hessian object)', stage_a)
log = pyscf.lib.logger.new_logger(hess, 6)

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

def repaired(h):
    """Diagonal-block translational sum-rule repair used by the production inputs."""
    fixed = h.copy()
    for p in range(h.shape[0]):
        fixed[p, p] = -(h[p].sum(axis=0) - h[p, p])
    return 0.5 * (fixed + fixed.transpose(1, 0, 3, 2))


def all_freqs(hessobj, h):
    res = hessobj.harmonic_analysis(mol, h, exclude_trans=False, exclude_rot=False,
                                    imaginary_freq=False, intensity=False)
    return numpy.sort(numpy.asarray(res['freq_wavenumber'], dtype=float))


def report(tag, hessobj, h):
    viol = numpy.abs(h.sum(axis=1)).reshape(mol.natm, -1).max(axis=1)
    print(f'{tag}: translational sum rule per atom (Eh/Bohr^2): '
          + ', '.join(f'{mol.atom_pure_symbol(i)}{i} {v:.2e}' for i, v in enumerate(viol)))
    print(f'{tag}: all 3N frequencies, raw Hessian (cm^-1):      '
          + ', '.join(f'{f:.1f}' for f in all_freqs(hessobj, h)))
    print(f'{tag}: all 3N frequencies, repaired Hessian (cm^-1): '
          + ', '.join(f'{f:.1f}' for f in all_freqs(hessobj, repaired(h))), flush=True)
    return viol


report('level 5', hess, de)
numpy.savez('hcnhbr_avtz_probe_hessian.npz', hessian_raw=de, hessian_repaired=repaired(de),
            coords=mol.atom_coords(), energy=e_scf)
stamp('STAGES A-F OK')

# ---------------------------------------------------------------------------
# Stage G: is the bromine sum-rule defect a grid-quadrature error?
#
# At level 5 bromine gets 120 radial x 770 angular points, pruned near the
# core by nwchem_prune, where the second derivatives of the core density are
# largest.  Each variant below re-runs SCF + full Hessian (the Hessian must be
# consistent with the SCF grid) and prints the same three lines as above.  If
# the Br violation falls steeply with the radial count, the defect is
# quadrature error and a denser Br grid is the fix; if it stays near 0.7, a
# term is missing in the XC Hessian and the repair stays the right tool.
# ---------------------------------------------------------------------------
GRID_VARIANTS = [
    ('level 5, no pruning',                lambda g: setattr(g, 'prune', None)),
    ('level 5, Br 240x770',                lambda g: setattr(g, 'atom_grid', {'Br': (240, 770)})),
    ('level 5, Br 240x770, no pruning',    lambda g: (setattr(g, 'atom_grid', {'Br': (240, 770)}),
                                                      setattr(g, 'prune', None))),
    ('level 5, Br 400x974, no pruning',    lambda g: (setattr(g, 'atom_grid', {'Br': (400, 974)}),
                                                      setattr(g, 'prune', None))),
]


def hessian_on_grid(tag, setup):
    m = neo.CDFT(mol, xc=XC, epc=None)
    m.verbose = 3
    m.conv_tol = 1e-11
    m.conv_tol_grad = 1e-6
    m.max_cycle = 200
    m.conv_tol_cpscf = 1e-9
    g = m.components['e'].grids
    g.level = GRID_LEVEL
    setup(g)
    e = m.kernel()
    print(f'{tag}: E = {e:.10f} Eh (level-5 default {e_scf:.10f}, difference {(e - e_scf) * 1e6:.3f} uEh), '
          f'converged = {m.converged}, grid points {g.weights.size}', flush=True)
    h = m.Hessian()
    h.max_cycle = HESS_MAX_CYCLE
    h.verbose = 4
    H = h.kernel()
    report(tag, h, H)
    return H


for tag, setup in GRID_VARIANTS:
    run_stage(f'G ({tag})', lambda tag=tag, setup=setup: hessian_on_grid(tag, setup))

stamp('ALL STAGES OK')
print_summary()
