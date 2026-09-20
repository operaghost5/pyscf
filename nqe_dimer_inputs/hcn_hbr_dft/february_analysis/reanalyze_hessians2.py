"""Re-run the harmonic analysis on the Hessians printed in PySCF CNEO logs (generic atoms).

Pure Python. For every log with a complete Hessian block:
  * parse the (natm,natm,3,3) Cartesian Hessian (Eh/Bohr^2), the geometry (Angstrom)
    and the atom symbols
  * mass-weight, diagonalise the full 3N matrix (no projection)
  * project out translations + rotations, dropping the on-axis rotation of a
    linear molecule (relative-norm test), and report the 3N-5 / 3N-6 vibrations
Frequencies in cm^-1; imaginary modes are printed with a trailing 'i'.
"""
import glob
import math
import os
import re
import sys

LOGDIR = sys.argv[1] if len(sys.argv) > 1 else "."

# PySCF isotope-averaged masses (pyscf.data.elements.MASSES)
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
        "Cl": 35.45, "Br": 79.904}
HARTREE2J = 4.3597447222071e-18
AMU = 1.66053906660e-27
BOHR_SI = 5.29177210903e-11
C_SI = 299792458.0
ANG2BOHR = 1e-10 / BOHR_SI
AU2CM = math.sqrt(HARTREE2J / (AMU * BOHR_SI ** 2)) / (2 * math.pi) / C_SI * 1e-2

BASIS = {"bas_two": "aVDZ", "bas_three": "aVTZ", "bas_four": "aVQZ", "bas_five": "aV5Z"}
FUNC = {"func_one": "PBE", "func_two": "PW91", "func_three": "BP86",
        "func_four": "BLYP", "func_five": "B97"}
REAL = re.compile(r"[-+]?\d+\.\d+(?:e[-+]?\d+)?|[-+]?\d+e[-+]?\d+")


def jacobi_eigvals(A, tol=1e-24, maxsweeps=200):
    n = len(A)
    A = [row[:] for row in A]
    for _ in range(maxsweeps):
        off = sum(A[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(A[p][q]) < 1e-300:
                    continue
                theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q])
                t = (1.0 if theta >= 0 else -1.0) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = A[k][p], A[k][q]
                    A[k][p] = c * akp - s * akq
                    A[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = A[p][k], A[q][k]
                    A[p][k] = c * apk - s * aqk
                    A[q][k] = s * apk + c * aqk
    return sorted(A[i][i] for i in range(n))


def matmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    return [[sum(A[i][l] * B[l][j] for l in range(k)) for j in range(m)] for i in range(n)]


def gram_schmidt(vectors, rel_drop_tol=1e-3):
    scale = max(math.sqrt(sum(x * x for x in v)) for v in vectors)
    basis, dropped = [], []
    for idx, v in enumerate(vectors):
        w = v[:]
        for b in basis:
            d = sum(x * y for x, y in zip(w, b))
            w = [x - d * y for x, y in zip(w, b)]
        nrm = math.sqrt(sum(x * x for x in w))
        if nrm < rel_drop_tol * scale:
            dropped.append(idx)
            continue
        basis.append([x / nrm for x in w])
    return basis, dropped


def to_freq(lam):
    return math.sqrt(abs(lam)) * AU2CM * (1 if lam >= 0 else -1)


def fmt(fs):
    return ", ".join((f"{-x:.1f}i" if x < 0 else f"{x:.1f}") for x in fs)


def parse_symbols(text, natm_hint=None):
    m = re.search(r"\[\['([A-Z][a-z]?)', array", text)
    if m:
        syms = re.findall(r"\['([A-Z][a-z]?)', array\(", text)
        if syms:
            return syms
    # fall back: last per-cycle geometry print
    geos = re.findall(r"^cycle \d+: E = .*\n((?:[A-Z][a-z]?\s+\S+\s+\S+\s+\S+\n)+)", text, re.M)
    if geos:
        return [l.split()[0] for l in geos[-1].strip().splitlines()]
    return None


def parse_hessian(text, natm):
    i0 = text.find("BEGIN HESSIAN CALCULATION")
    if i0 < 0:
        return None
    a = text.find("[[[[", i0)
    b = text.find("]]]]", a)
    if a < 0 or b < 0:
        return None
    vals = [float(x) for x in REAL.findall(text[a:b + 4])]
    n3 = 3 * natm
    if len(vals) != n3 * n3:
        return None
    H = [[0.0] * n3 for _ in range(n3)]
    k = 0
    for p in range(natm):
        for q in range(natm):
            for x in range(3):
                for y in range(3):
                    H[3 * p + x][3 * q + y] = vals[k]
                    k += 1
    return H


def parse_geom(text, natm):
    for pat in (r"FINAL OPTIMIZED GEOMETRY: \n(\[\[.*?\]\])",
                r"BEGIN OPTIMIZED ATOMIC COORDINATES:\n-+\n(\[\[.*?\]\])"):
        m = re.search(pat, text, re.S)
        if m:
            nums = [float(x) for x in REAL.findall(m.group(1))]
            if len(nums) == 3 * natm:
                return [nums[i:i + 3] for i in range(0, 3 * natm, 3)]
    return None


def analyse(H, xyz_ang, symbols):
    natm = len(symbols)
    n3 = 3 * natm
    mass = [MASS[s] for s in symbols]
    asym = max(abs(H[i][j] - H[j][i]) for i in range(n3) for j in range(n3))
    Hs = [[0.5 * (H[i][j] + H[j][i]) for j in range(n3)] for i in range(n3)]
    sm = [math.sqrt(mass[i // 3]) for i in range(n3)]
    Hm = [[Hs[i][j] / (sm[i] * sm[j]) for j in range(n3)] for i in range(n3)]
    full = [to_freq(l) for l in jacobi_eigvals(Hm)]

    xyz = [[c * ANG2BOHR for c in r] for r in xyz_ang]
    mtot = sum(mass)
    com = [sum(mass[i] * xyz[i][k] for i in range(natm)) / mtot for k in range(3)]
    r = [[xyz[i][k] - com[k] for k in range(3)] for i in range(natm)]
    TR = []
    for ax in range(3):
        v = [0.0] * n3
        for i in range(natm):
            v[3 * i + ax] = math.sqrt(mass[i])
        TR.append(v)
    for ax in range(3):
        e = [0.0, 0.0, 0.0]
        e[ax] = 1.0
        v = [0.0] * n3
        for i in range(natm):
            cx = e[1] * r[i][2] - e[2] * r[i][1]
            cy = e[2] * r[i][0] - e[0] * r[i][2]
            cz = e[0] * r[i][1] - e[1] * r[i][0]
            v[3 * i:3 * i + 3] = [math.sqrt(mass[i]) * cx, math.sqrt(mass[i]) * cy, math.sqrt(mass[i]) * cz]
        TR.append(v)
    basis, dropped = gram_schmidt(TR)
    n_ext = len(basis)
    P = [[(1.0 if i == j else 0.0) - sum(b[i] * b[j] for b in basis) for j in range(n3)] for i in range(n3)]
    Hp = matmul(matmul(P, Hm), P)
    proj = [to_freq(l) for l in jacobi_eigvals(Hp)]
    proj_sorted = sorted(proj, key=abs)
    zeros, vib = proj_sorted[:n_ext], sorted(proj_sorted[n_ext:])
    # off-axis distance of any atom from the line through the two heaviest-separated atoms (0 and -1)
    a, b = xyz_ang[0], xyz_ang[-1]
    d = [b[k] - a[k] for k in range(3)]
    dn = math.sqrt(sum(x * x for x in d))
    d = [x / dn for x in d]
    offaxis = 0.0
    for p in xyz_ang:
        v = [p[k] - a[k] for k in range(3)]
        t = sum(v[k] * d[k] for k in range(3))
        offaxis = max(offaxis, math.sqrt(max(0.0, sum(x * x for x in v) - t * t)))
    return dict(full=full, vib=vib, zeros=zeros, n_ext=n_ext, linear=(n_ext == 5),
                asym=asym, offaxis=offaxis)


logs = sorted(glob.glob(os.path.join(LOGDIR, "*.log")))
for p in logs:
    t = open(p, errors="replace").read()
    if "END HESSIAN CALCULATION" not in t:
        continue
    symbols = parse_symbols(t)
    if not symbols:
        print(f"{p}: cannot determine atom symbols")
        continue
    natm = len(symbols)
    H = parse_hessian(t, natm)
    g = parse_geom(t, natm)
    if H is None or g is None:
        print(f"{p}: Hessian or geometry not parseable")
        continue
    parts = os.path.basename(p).split(".")
    logged = re.search(r"BEGIN HARMONIC ANALYSIS.*?\n-+\n(.*?)\n-+\nEND HARMONIC", t, re.S)
    res = analyse(H, g, symbols)
    print("=" * 96)
    print(f"{parts[4]:10s} {FUNC.get(parts[3], parts[3]):5s} {BASIS.get(parts[1], parts[1]):5s} {parts[5]}   atoms {' '.join(symbols)}")
    print(f"  off-axis {res['offaxis']:.1e} Ang | max|H-H^T| {res['asym']:.1e} | "
          f"external modes projected: {res['n_ext']} ({'linear' if res['linear'] else 'non-linear'})")
    if logged:
        print("  PySCF logged      : " + re.sub(r"\s+", " ", logged.group(1)).strip()[:400])
    print("  full 3N (no proj) : " + fmt(res['full']))
    nimag = sum(1 for f in res['vib'] if f < 0)
    print(f"  VIBRATIONAL ({len(res['vib'])})  : " + fmt(res['vib']) + f"    -> imaginary modes: {nimag}")
