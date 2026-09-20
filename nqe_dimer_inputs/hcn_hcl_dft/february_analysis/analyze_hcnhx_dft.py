"""HCN...HX (X = F, Cl, Br): classification of the February 2026 conventional-DFT runs (aug-cc-pV5Z excluded).

Reuses the parsing / projection functions of reanalyze_hessians2.py (copied alongside) and
adds the translational sum rule per atom, the diagonal-block repaired spectrum, the geometry
descriptors (hydrogen-bonded linear, hydrogen-bonded bent or halogen-bonded linear, read from the
distances), the optimizer status and, when a directory of CNEO-DFT logs is given as the second
argument, the CNEO - DFT frequency and geometry shifts for the runs that have a Hessian in both.

    python analyze_hcnhx_dft.py <DFT log dir> [<CNEO-DFT log dir>]

Atom order in every log: H0 C1 N2 H3 X4.  The halogen X is read from the logs.  Pure Python.
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
    records, dups = {}, []
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
                   kry_max=(max(kry) + 1) if kry else None, file=os.path.basename(p),
                   last_line=t.strip().splitlines()[-1] if t.strip() else "")
        m = re.search(r"FINAL SINGLE POINT ENERGY:\s+(\S+)", t)
        rec["e_sp"] = float(m.group(1)) if m else None
        if rec["hess_end"]:
            symbols = parse_symbols(t)
            H = parse_hessian(t, len(symbols)) if symbols else None
            g = parse_geom(t, len(symbols)) if symbols else None
            if symbols is None or H is None or g is None:
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
        key = (start, rec["basis"], rec["xc"])
        old = records.get(key)
        if old is None or (rec["hess_end"] and not old["hess_end"]) or \
           (rec["hess_end"] == old["hess_end"] and date_key(date) > date_key(old["date"])):
            if old is not None:
                dups.append(f"{key}: using {rec['file']} ({rec['ncyc']} steps) over {old['file']} ({old['ncyc']} steps)")
            records[key] = rec
        elif old is not None:
            dups.append(f"{key}: using {old['file']} ({old['ncyc']} steps) over {rec['file']} ({rec['ncyc']} steps)")
    return records, dups


def date_key(d):
    """'07-07-2026' -> (2026, 7, 7) so that later reruns sort after earlier ones."""
    try:
        dd, mm, yy = d.split("-")
        return (int(yy), int(mm), int(dd))
    except ValueError:
        return (0, 0, 0)


def descriptors(x, X):
    """H0 C1 N2 H3 X4.  Distances in Angstrom, angles in degrees; structure type from the distances."""
    d = dict(rCH=dist(x[0], x[1]), rCN=dist(x[1], x[2]), rHX=dist(x[3], x[4]),
             rN_H=dist(x[2], x[3]), rN_X=dist(x[2], x[4]), rH0_X=dist(x[0], x[4]), rC_X=dist(x[1], x[4]))
    contacts = {"HX": d["rN_H"], "XB": d["rN_X"], "CH": d["rH0_X"]}
    kind = min(contacts, key=contacts.get)
    if kind == "HX":                                        # N...H-X, HX the donor
        d["type"] = f"N...H-{X}"
        d["a1"], d["a2"] = angle(x[4], x[3], x[2]), angle(x[3], x[2], x[1])      # X-H...N, H...N-C
    elif kind == "XB":                                      # N...X-H halogen bond
        d["type"] = f"N...{X}-H"
        d["a1"], d["a2"] = angle(x[2], x[4], x[3]), angle(x[4], x[2], x[1])      # N...X-H, X...N-C
    else:                                                   # C-H...X, the C-H the donor
        d["type"] = f"C-H...{X}"
        d["a1"], d["a2"] = angle(x[1], x[0], x[4]), angle(x[0], x[4], x[3])      # C-H...X, H...X-H
    d["closest"] = min(contacts.values())
    d["bound"] = d["closest"] < 4.0
    return d


records, dups = read_logs(LOGDIR, "dft")
cneo, cdups = read_logs(CNEODIR, "cneodft") if CNEODIR else ({}, [])
X = None
for r in records.values():
    if r.get("symbols"):
        X = r["symbols"][4]
        break
if X is None:
    sys.exit("no DFT log with a Hessian found")

out = []
w = out.append
w(f"HCN...H{X}, conventional DFT (all nuclei classical), February 2026 runs: classification (aug-cc-pV5Z excluded)")
w(f"Masses: PySCF isotope-averaged (H 1.008, C 12.011, N 14.007, {X} {MASS[X]}).  Atom order in every log: H0 C1 N2 H3 {X}4.")
w("Original settings (from the .py inputs): dft.KS, grids.level = 9 (200 radial shells x 1454 angular points on every atom),")
w("SCF conv_tol 1e-9 (default), geomeTRIC with the custom criteria (dE 1e-9, grms 3e-6, gmax 4.5e-4, drms 1.2e-3, dmax 1.8e-3;")
w("maxiter 1500), analytic RKS Hessian with the default CPHF (50 Krylov iterations), thermo.harmonic_analysis with default projection.")
if dups:
    w("Duplicate logs for one case (the later date wins when both have a Hessian):")
    for d in dups:
        w("  " + d)
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
            if r.get("note"):
                st += "; " + r["note"]
            w(f"{start:<11}{xc:<6}{bas:<13}{r['ncyc']:>9d} {str(r['kry_max'] or '-'):>8}  {st}")
w("")
w("PART 1: geometry, energy and classification of the runs that produced a Hessian (raw Hessian, correct 3N-5/3N-6 projection)")
w(f"  type = the shortest intermolecular contact: N...H-{X} (H{X} the donor, linear), C-H...{X} (the C-H the donor, bent) or")
w(f"  N...{X}-H (halogen bond, linear); a1/a2 = the X-H...Y and the acceptor angle of that contact (F-H...N & H...N-C, C-H...X & H...X-H,")
w(f"  or N...X-H & X...N-C); off-axis = largest distance of any atom from the H0-{X}4 line; '{X} viol' = {X} sum-rule violation (Eh/Bohr^2).")
w(f"{'start':<11}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'C-H':>6} {'C-N':>6} {'H-'+X:>6} {'N..H':>6} {'N..'+X:>6} {'H0..'+X:>6} {'type':>9} {'a1':>6} {'a2':>6} "
  f"{'offaxis':>8} {X+' viol':>8}  vibrations (cm^-1)                                                           imag  class")
class_summary = []
for start in STARTS:
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None or not r.get("hess_end"):
                continue
            d = descriptors(r["xyz"], X)
            vib = r["res"]["vib"]
            nimag = sum(1 for f in vib if f < 0)
            cls = "MINIMUM" if nimag == 0 else f"SADDLE (order {nimag})"
            flag = "" if d["bound"] else "  UNBOUND (closest contact > 4 A)"
            w(f"{start:<11}{xc:<6}{bas:<13}{r['e_sp']:>16.9f} {d['rCH']:6.3f} {d['rCN']:6.3f} {d['rHX']:6.3f} {d['rN_H']:6.3f} {d['rN_X']:6.3f} {d['rH0_X']:6.3f} "
              f"{d['type']:>9} {d['a1']:6.1f} {d['a2']:6.1f} {r['res']['offaxis']:8.1e} {r['viol'][4]:8.1e}  {fmt(vib):<75s} {nimag:>3d}  {cls}{flag}")
            class_summary.append((start, bas, xc, cls, r["e_sp"], d["type"], d["bound"]))
w("")
w("PART 2: external modes, the sum rule and the repaired spectrum; PySCF's own harmonic analysis as logged")
w("  'full 3N' is the unprojected spectrum (5 external modes for a linear complex, 6 for the bent one, all should be ~0);")
w("  'repaired' resets each diagonal block to minus the sum of its off-diagonal blocks.  The logged PySCF list is what the")
w("  pre-2.14 fork printed; for a linear complex it has 9 entries where 10 (3N-5) are expected when the rotor test misfired.")
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
w("  energies relative to the global_min run of the same functional/basis (kcal/mol); all nuclei classical, so directly comparable:")
w(f"    {'xc':<6}{'basis':<13}{'local_min':>10}{'xlocal_min':>11}   structure types (global / local / xlocal)")
for bas in BAS_ORDER:
    for xc in FUNC_ORDER:
        rows = {c[0]: c for c in class_summary if c[1] == bas and c[2] == xc}
        if "global_min" not in rows:
            continue
        eg = rows["global_min"][4]
        el = f"{(rows['local_min'][4] - eg) * K:+10.2f}" if "local_min" in rows else f"{'-':>10}"
        ex = f"{(rows['xlocal_min'][4] - eg) * K:+11.2f}" if "xlocal_min" in rows else f"{'-':>11}"
        types = " / ".join(rows[s][5] + ("" if rows[s][6] else " (unbound)") if s in rows else "-" for s in STARTS)
        w(f"    {xc:<6}{bas:<13}{el}{ex}   {types}")

if cneo:
    w("")
    w("PART 4: CNEO-DFT minus DFT, same functional, basis and start, for the runs with a Hessian in both sets")
    w("  Frequencies from the repaired Hessians, same masses (H 1.008) on both sides so the shift is the method alone;")
    w("  modes matched by rank after sorting.  Geometry shifts in Angstrom (CNEO - DFT).")
    if cdups:
        w("  duplicate CNEO logs: " + "; ".join(cdups))
    for start in STARTS:
        for bas in BAS_ORDER:
            for xc in FUNC_ORDER:
                r = records.get((start, bas, xc))
                c = cneo.get((start, bas, xc))
                if r is None or c is None or not r.get("hess_end") or not c.get("hess_end"):
                    continue
                dr, dc = descriptors(r["xyz"], X), descriptors(c["xyz"], X)
                w("=" * 110)
                w(f"{start} {xc} {bas}   DFT structure {dr['type']}{'' if dr['bound'] else ' (unbound)'}, "
                  f"CNEO structure {dc['type']}{'' if dc['bound'] else ' (unbound)'}")
                w(f"  DFT  vibrations: {fmt(r['resR']['vib'])}")
                w(f"  CNEO vibrations: {fmt(c['resR']['vib'])}")
                if len(r["resR"]["vib"]) == len(c["resR"]["vib"]):
                    w("  CNEO - DFT     : " + ", ".join(f"{b - a:+.1f}" for a, b in zip(r["resR"]["vib"], c["resR"]["vib"])))
                else:
                    w("  (different numbers of vibrations: one structure is linear and the other is not; no mode-by-mode shift)")
                w(f"  CNEO - DFT geometry: d(C-H) {dc['rCH'] - dr['rCH']:+.4f}  d(C-N) {dc['rCN'] - dr['rCN']:+.4f}  d(H-{X}) {dc['rHX'] - dr['rHX']:+.4f}  "
                  f"d(N..H) {dc['rN_H'] - dr['rN_H']:+.4f}  d(N..{X}) {dc['rN_X'] - dr['rN_X']:+.4f}  d(H0..{X}) {dc['rH0_X'] - dr['rH0_X']:+.4f}")
                w(f"  CNEO - DFT energy (Eh): {c['e_sp'] - r['e_sp']:+.6f}   (includes the protons' zero-point-like energy; not a binding-energy shift)")

text = "\n".join(out) + "\n"
open(os.path.join(HERE, f"hcn_h{X.lower()}_dft_classification.txt"), "w").write(text)
print(text)
