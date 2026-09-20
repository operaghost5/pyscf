"""HCN...HF: classification of the February 2026 conventional-DFT runs (aug-cc-pV5Z excluded).

Reuses the parsing / projection functions of reanalyze_hessians2.py (copied alongside) and
adds the translational sum rule per atom, the diagonal-block repaired spectrum, the geometry
descriptors, the optimizer status and, when a directory of CNEO-DFT logs is given as the
second argument, the CNEO - DFT frequency shifts for the runs that have a Hessian in both.

    python analyze_hcnhf_dft.py <DFT log dir> [<CNEO-DFT log dir>]

Pure Python.
"""
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, "reanalyze_hessians2.py")).read()
exec(SRC.split("\nlogs = sorted(")[0])          # functions and constants only, no main loop

LOGDIR = sys.argv[1]
CNEODIR = sys.argv[2] if len(sys.argv) > 2 else None
BASNAME = {"bas_two": "aug-cc-pVDZ", "bas_three": "aug-cc-pVTZ", "bas_four": "aug-cc-pVQZ"}
FUNC_ORDER = ["PBE", "PW91", "BP86", "BLYP", "B97"]
BAS_ORDER = ["aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]
STARTS = ("global_min", "local_min", "xlocal_min")
K = 627.5095


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def angle(a, b, c):
    """angle a-b-c in degrees"""
    u = [a[k] - b[k] for k in range(3)]
    v = [c[k] - b[k] for k in range(3)]
    cs = sum(u[k] * v[k] for k in range(3)) / (dist(a, b) * dist(c, b))
    return math.degrees(math.acos(max(-1.0, min(1.0, cs))))


def sum_rule(H, natm):
    viol = []
    for p in range(natm):
        worst = 0.0
        for x in range(3):
            for y in range(3):
                s = sum(H[3 * p + x][3 * q + y] for q in range(natm))
                worst = max(worst, abs(s))
        viol.append(worst)
    return viol


def repair(H, natm):
    n3 = 3 * natm
    R = [row[:] for row in H]
    for p in range(natm):
        for x in range(3):
            for y in range(3):
                off = sum(H[3 * p + x][3 * q + y] for q in range(natm) if q != p)
                R[3 * p + x][3 * p + y] = -off
    return [[0.5 * (R[i][j] + R[j][i]) for j in range(n3)] for i in range(n3)]


def read_logs(logdir, method_token):
    records = {}
    for p in sorted(glob.glob(os.path.join(logdir, "*.log"))):
        parts = os.path.basename(p).split(".")
        if len(parts) < 6 or parts[2] != method_token:
            continue
        bas, func, start, date = parts[1], parts[3], parts[4], parts[5]
        if bas not in BASNAME:
            continue
        t = open(p, errors="replace").read()
        cycles = re.findall(r"^cycle (\d+): E = (\S+)\s+dE = (\S+)\s+norm\(grad\) = (\S+)", t, re.M)
        kry = [int(x) for x in re.findall(r"krylov cycle (\d+)", t[t.find("BEGIN HESSIAN CALCULATION"):])] \
            if "BEGIN HESSIAN CALCULATION" in t else []
        rec = dict(start=start, xc=FUNC[func], basis=BASNAME[bas], date=date, ncyc=len(cycles),
                   lines=t.count("\n"), scf_started="BEGIN INITIAL SCF RUN" in t,
                   opt_end="END GEOMETRIC GEOMETRY OPTIMIZATION" in t,
                   capped="Geometry optimization is not converged" in t,
                   hess_end="END HESSIAN CALCULATION" in t,
                   scf_unconv=len(re.findall(r"SCF not converged", t)),
                   kry_max=(max(kry) + 1) if kry else None,
                   last_line=t.strip().splitlines()[-1] if t.strip() else "")
        m = re.search(r"FINAL SINGLE POINT ENERGY:\s+(\S+)", t)
        rec["e_sp"] = float(m.group(1)) if m else None
        if rec["hess_end"]:
            symbols = parse_symbols(t)
            H = parse_hessian(t, len(symbols)) if symbols else None
            g = parse_geom(t, len(symbols)) if symbols else None
            if symbols is None or H is None or g is None:
                # truncated log (e.g. a rerun that died while printing the summary): keep any complete record
                rec["hess_end"] = False
                rec["note"] = "Hessian printed but the summary geometry is missing; log truncated"
        if rec["hess_end"]:
            natm = len(symbols)
            res = analyse(H, g, symbols)
            viol = sum_rule(H, natm)
            resR = analyse(repair(H, natm), g, symbols)
            logged = re.search(r"BEGIN HARMONIC ANALYSIS.*?\n-+\n(.*?)\n-+\nEND HARMONIC", t, re.S)
            rec.update(symbols=symbols, xyz=g, res=res, resR=resR, viol=viol,
                       logged=[float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:e[-+]?\d+)?", logged.group(1))]
                       if logged else None)
        # a rerun of the same case on a later date replaces the earlier one only if it got further
        key = (start, rec["basis"], rec["xc"])
        old = records.get(key)
        if old is None or (rec["hess_end"] and not old["hess_end"]) or \
           (rec["hess_end"] == old["hess_end"] and date > old["date"]):
            records[key] = rec
    return records


def descriptors(x):
    """H0 C1 N2 H3 F4.  Returns dict of distances (Angstrom) and angles (deg)."""
    d = dict(rCH=dist(x[0], x[1]), rCN=dist(x[1], x[2]), rHF=dist(x[3], x[4]),
             rN_H=dist(x[2], x[3]), rN_F=dist(x[2], x[4]),
             rH0_F=dist(x[0], x[4]), rC_F=dist(x[1], x[4]))
    if d["rN_H"] < d["rH0_F"]:
        d["donor"] = "HF"                                  # N...H-F
        d["a_XH_Y"] = angle(x[4], x[3], x[2])              # F-H...N
        d["a_HYZ"] = angle(x[3], x[2], x[1])               # H...N-C
    else:
        d["donor"] = "CH"                                  # C-H...F
        d["a_XH_Y"] = angle(x[1], x[0], x[4])              # C-H...F
        d["a_HYZ"] = angle(x[0], x[4], x[3])               # H...F-H
    return d


records = read_logs(LOGDIR, "dft")
cneo = read_logs(CNEODIR, "cneodft") if CNEODIR else {}

out = []
w = out.append
w("HCN...HF, conventional DFT (all nuclei classical), February 2026 runs: classification (aug-cc-pV5Z excluded)")
w("Masses: PySCF isotope-averaged (H 1.008, C 12.011, N 14.007, F 18.998).  Atom order in every log: H0 C1 N2 H3 F4.")
w("Original settings (from the .py inputs): dft.KS, grids.level = 9 (200 radial shells x 1454 angular points on every atom),")
w("SCF conv_tol 1e-9 (default), geomeTRIC with the custom criteria (dE 1e-9, grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3;")
w("maxiter 1500), analytic RKS Hessian with the default CPHF (50 Krylov iterations), thermo.harmonic_analysis with default projection.")
w("")
w("PART 0: run status")
w(f"{'start':<11}{'xc':<6}{'basis':<13}{'opt steps':>9} {'CPHF it':>8}  status")
for start in STARTS:
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None:
                w(f"{start:<11}{xc:<6}{bas:<13}{'-':>9} {'-':>8}  no log")
                continue
            if not r["scf_started"]:
                st = f"DIED BEFORE THE FIRST SCF ({r['lines']} lines; last line: '{r['last_line']}')"
            elif not r["opt_end"]:
                st = "optimization did not finish"
            elif r["capped"]:
                st = "optimizer hit its step cap; Hessian run anyway" if r["hess_end"] else "optimizer hit its step cap; no Hessian"
            else:
                st = "optimizer converged; Hessian done" if r["hess_end"] else "optimizer converged; no Hessian"
            if r["scf_unconv"]:
                st += f"; {r['scf_unconv']} 'SCF not converged' message(s)"
            w(f"{start:<11}{xc:<6}{bas:<13}{r['ncyc']:>9d} {str(r['kry_max'] or '-'):>8}  {st}")
w("")
w("PART 1: geometry, energy and classification of the runs that produced a Hessian (raw Hessian, correct 3N-5/3N-6 projection)")
w("  donor = HF when N...H(F) is shorter than H(C)...F (linear N...H-F complex), CH otherwise (bent C-H...F-H complex);")
w("  XH..Y = the X-H...Y angle of that hydrogen bond, H..YZ = the angle at the acceptor (H...N-C or H...F-H); off-axis = largest")
w("  distance of any atom from the H0-F4 line; 'F viol' = fluorine sum-rule violation (Eh/Bohr^2).")
w(f"{'start':<11}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'C-H':>6} {'C-N':>6} {'H-F':>6} {'N..H':>6} {'H0..F':>6} {'N..F':>6} {'donor':>5} {'XH..Y':>6} {'H..YZ':>6} "
  f"{'offaxis':>8} {'F viol':>8}  vibrations (cm^-1)                                                          imag  class")
class_summary = []
for start in STARTS:
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None or not r.get("hess_end"):
                continue
            d = descriptors(r["xyz"])
            vib = r["res"]["vib"]
            nimag = sum(1 for f in vib if f < 0)
            cls = "MINIMUM" if nimag == 0 else f"SADDLE (order {nimag})"
            w(f"{start:<11}{xc:<6}{bas:<13}{r['e_sp']:>16.9f} {d['rCH']:6.3f} {d['rCN']:6.3f} {d['rHF']:6.3f} {d['rN_H']:6.3f} {d['rH0_F']:6.3f} {d['rN_F']:6.3f} "
              f"{d['donor']:>5} {d['a_XH_Y']:6.1f} {d['a_HYZ']:6.1f} {r['res']['offaxis']:8.1e} {r['viol'][4]:8.1e}  {fmt(vib):<75s} {nimag:>3d}  {cls}")
            class_summary.append((start, bas, xc, cls, r["e_sp"], d["donor"]))
w("")
w("PART 2: external modes, the sum rule and the repaired spectrum; PySCF's own harmonic analysis as logged")
w("  'full 3N' is the unprojected spectrum (5 external modes for the linear complex, 6 for the bent one, all should be ~0);")
w("  'repaired' resets each diagonal block to minus the sum of its off-diagonal blocks.  The logged PySCF list is what the")
w("  pre-2.14 fork printed; for the linear complex it has 9 entries where 10 (3N-5) are expected when the rotor test misfired.")
for start in STARTS:
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None or not r.get("hess_end"):
                continue
            v = r["viol"]
            w("=" * 110)
            w(f"{start} {xc} {bas}   {'linear' if r['res']['linear'] else 'non-linear'} ({r['res']['n_ext']} external modes); "
              "sum rule per atom (Eh/Bohr^2): " + ", ".join(f"{s}{i} {vv:.1e}" for i, (s, vv) in enumerate(zip(r["symbols"], v))))
            lg = r["logged"] or []
            w(f"  PySCF logged ({len(lg)} values)   : " + fmt(lg))
            w("  re-diagonalized, projected: " + fmt(r["res"]["vib"]))
            w("  full 3N, raw              : " + fmt(r["res"]["full"]))
            w("  full 3N, repaired         : " + fmt(r["resR"]["full"]))
            w("  repaired, projected       : " + fmt(r["resR"]["vib"]))
            dmax = max(abs(a - b) for a, b in zip(r["res"]["vib"], r["resR"]["vib"]))
            w(f"  max |raw - repaired| over the vibrations: {dmax:.2f} cm^-1")
w("")
w("PART 3: classification on the RAW and on the REPAIRED Hessian")
w(f"{'start':<11}{'xc':<6}{'basis':<13} {'raw class':<18}{'repaired class':<16}  repaired vibrations (cm^-1)")
for start in STARTS:
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None or not r.get("hess_end"):
                continue
            vr, vR = r["res"]["vib"], r["resR"]["vib"]
            nr, nR = sum(1 for f in vr if f < 0), sum(1 for f in vR if f < 0)
            cr = "MINIMUM" if nr == 0 else f"SADDLE (order {nr})"
            cR = "MINIMUM" if nR == 0 else f"SADDLE (order {nR})"
            w(f"{start:<11}{xc:<6}{bas:<13} {cr:<18}{cR:<16}  {fmt(vR)}")
w("")
w("SUMMARY")
n_min = sum(1 for c in class_summary if c[3] == "MINIMUM")
nR_min = sum(1 for k, r in records.items() if r.get("hess_end") and not any(f < 0 for f in r["resR"]["vib"]))
w(f"  Hessians available: {len(class_summary)}.  Raw Hessian: {n_min} minima, {len(class_summary) - n_min} apparent saddle points."
  f"  Repaired Hessian: {nR_min} minima, {len(class_summary) - nR_min} saddle points.")
w("  energy of the local_min run relative to the global_min run of the same functional/basis (kcal/mol); all nuclei classical,")
w("  so the two isomers are directly comparable:")
for bas in BAS_ORDER:
    for xc in FUNC_ORDER:
        eg = [c[4] for c in class_summary if c[0] == "global_min" and c[1] == bas and c[2] == xc]
        el = [c[4] for c in class_summary if c[0] == "local_min" and c[1] == bas and c[2] == xc]
        if eg and el:
            w(f"    {xc:<6}{bas:<13}{(el[0] - eg[0]) * K:+8.2f}")

if cneo:
    w("")
    w("PART 4: CNEO-DFT minus DFT, same functional, basis and start, for the runs with a Hessian in both sets")
    w("  Frequencies from the repaired Hessians, same masses (H 1.008) on both sides so the shift is the method alone;")
    w("  modes matched by rank after sorting.  Geometry shifts in Angstrom (CNEO - DFT).")
    for start in STARTS:
        for bas in BAS_ORDER:
            for xc in FUNC_ORDER:
                r = records.get((start, bas, xc))
                c = cneo.get((start, bas, xc))
                if r is None or c is None or not r.get("hess_end") or not c.get("hess_end"):
                    continue
                dr, dc = descriptors(r["xyz"]), descriptors(c["xyz"])
                w("=" * 110)
                w(f"{start} {xc} {bas}")
                w(f"  DFT  vibrations: {fmt(r['resR']['vib'])}")
                w(f"  CNEO vibrations: {fmt(c['resR']['vib'])}")
                if len(r["resR"]["vib"]) == len(c["resR"]["vib"]):
                    w("  CNEO - DFT     : " + ", ".join(f"{b - a:+.1f}" for a, b in zip(r["resR"]["vib"], c["resR"]["vib"])))
                w(f"  CNEO - DFT geometry: d(C-H) {dc['rCH'] - dr['rCH']:+.4f}  d(C-N) {dc['rCN'] - dr['rCN']:+.4f}  d(H-F) {dc['rHF'] - dr['rHF']:+.4f}  "
                  f"d(N..H) {dc['rN_H'] - dr['rN_H']:+.4f}  d(N..F) {dc['rN_F'] - dr['rN_F']:+.4f}")
                w(f"  CNEO - DFT energy (Eh): {c['e_sp'] - r['e_sp']:+.6f}   (includes the protons' zero-point-like energy; not a binding-energy shift)")

text = "\n".join(out) + "\n"
open(os.path.join(HERE, "hcn_hf_dft_classification.txt"), "w").write(text)
print(text)
