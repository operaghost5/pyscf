"""(HX)2 dimers (X = F, Cl, Br): classification of the April 2026 conventional-DFT runs (aug-cc-pV5Z excluded).

Reuses the parsing / projection functions of reanalyze_hessians2.py (copied alongside) and
adds the translational sum rule per atom, the diagonal-block repaired spectrum, the geometry
descriptors (covalent partners and the donor read from the distances), the optimizer status
and, when a directory of CNEO-DFT logs is given as the second argument, the CNEO - DFT
frequency and geometry shifts for the runs that have a Hessian in both, with the quantum
proton(s) of each CNEO run noted.

    python analyze_hxhx_dft.py <DFT log dir> [<CNEO-DFT log dir>]

Atom order in every log: H0 X1 H2 X3.  The halogen X is read from the logs.  Pure Python.
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
STARTS = ("global_min", "local_min")
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


def date_key(d):
    try:
        dd, mm, yy = d.split("-")
        return (int(yy), int(mm), int(dd))
    except ValueError:
        return (0, 0, 0)


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
                   qnuc=sorted(set(re.findall(r"\|g_(n\d+)\|=", t))),
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
                dups.append(f"{key}: using {rec['file']} over {old['file']}")
            records[key] = rec
        elif old is not None:
            dups.append(f"{key}: using {old['file']} over {rec['file']}")
    return records, dups


def descriptors(x, X):
    """H0 X1 H2 X3.  Covalent partners and the donor from the distances."""
    # covalent partner of each hydrogen: the nearer halogen
    partner = {h: (1 if dist(x[h], x[1]) < dist(x[h], x[3]) else 3) for h in (0, 2)}
    other = {1: 3, 3: 1}
    if partner[0] == partner[2]:                       # both H on one halogen: not an (HX)2 dimer
        return dict(valid=False)
    # hydrogen-bond candidates: each H to the other halogen
    hb = {h: dist(x[h], x[other[partner[h]]]) for h in (0, 2)}
    hd = min(hb, key=hb.get)                           # donor H
    ha = 2 if hd == 0 else 0                           # acceptor H
    Xd, Xa = partner[hd], partner[ha]
    d = dict(valid=True, hd=hd, ha=ha, Xd=Xd, Xa=Xa,
             rXHd=dist(x[hd], x[Xd]), rHX=hb[hd], rXHa=dist(x[ha], x[Xa]), rXX=dist(x[Xd], x[Xa]),
             rHXb=hb[ha],
             a_XHX=angle(x[Xd], x[hd], x[Xa]), a_HXH=angle(x[hd], x[Xa], x[ha]))
    d["asym"] = max(abs(d["rXHd"] - d["rXHa"]), abs(d["rHX"] - d["rHXb"]))
    d["label"] = f"H{hd}-{X}{Xd}...{X}{Xa}-H{ha}"
    d["bound"] = d["rHX"] < 4.0
    return d


records, dups = read_logs(LOGDIR, "dft")
cneo, cdups = read_logs(CNEODIR, "cneodft") if CNEODIR else ({}, [])
X = None
for r in records.values():
    if r.get("symbols"):
        X = r["symbols"][1]
        break
if X is None:
    sys.exit("no DFT log with a Hessian found")

out = []
w = out.append
w(f"(H{X})2, conventional DFT (all nuclei classical), April 2026 runs: classification (aug-cc-pV5Z excluded)")
w(f"Masses: PySCF isotope-averaged (H 1.008, {X} {MASS[X]}).  Atom order in every log: H0 {X}1 H2 {X}3; the covalent partners")
w("and the donor are read from the distances.")
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
w("PART 1: geometry, energy and classification of the runs that produced a Hessian (raw Hessian, 3N-6 projection)")
w(f"  bonding = donor H, its {X}, the acceptor {X} and its H, read from the distances; X-Hd = donor covalent bond, H..X = hydrogen bond,")
w(f"  X-Ha = acceptor covalent bond, X..X; angles X-H...X (linearity of the H bond) and H...X-H (tilt of the acceptor);")
w(f"  asym = C2h asymmetry, max(|X-Hd - X-Ha|, |H..X - H'..X'|); '{X} viol' = larger of the two {X} sum-rule violations (Eh/Bohr^2).")
w(f"{'start':<11}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'bonding':<15} {'X-Hd':>6} {'H..X':>6} {'X-Ha':>6} {'X..X':>6} {'XH..X':>6} {'H..XH':>6} {'asym':>6} "
  f"{X+' viol':>8}  vibrations (cm^-1), 3N-6 = 6                                      imag  class")
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
            fv = max(r["viol"][1], r["viol"][3])
            if not d["valid"]:
                w(f"{start:<11}{xc:<6}{bas:<13}{r['e_sp']:>16.9f} {'not an (HX)2 dimer':<15}")
                continue
            flag = "" if d["bound"] else "  UNBOUND (H..X > 4 A)"
            w(f"{start:<11}{xc:<6}{bas:<13}{r['e_sp']:>16.9f} {d['label']:<15} {d['rXHd']:6.3f} {d['rHX']:6.3f} {d['rXHa']:6.3f} {d['rXX']:6.3f} "
              f"{d['a_XHX']:6.1f} {d['a_HXH']:6.1f} {d['asym']:6.3f} {fv:8.1e}  {fmt(vib):<58s} {nimag:>3d}  {cls}{flag}")
            class_summary.append((start, bas, xc, cls, r["e_sp"], d["label"]))
w("")
w("PART 2: external modes, the sum rule and the repaired spectrum; PySCF's own harmonic analysis as logged")
w("  'full 3N' is the unprojected spectrum (6 external modes, all should be ~0); 'repaired' resets each diagonal block")
w("  to minus the sum of its off-diagonal blocks.")
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
n_local = sum(1 for k in records if k[0] == "local_min")
n_local_dead = sum(1 for k, r in records.items() if k[0] == "local_min" and not r["scf_started"])
w(f"  local_min logs: {n_local}, of which {n_local_dead} died before the first SCF (undefined 'DNE' starting geometry); nothing to classify there.")

if cneo:
    w("")
    w("PART 4: CNEO-DFT minus DFT, same functional and basis, for the runs with a Hessian in both sets")
    w("  Frequencies from the repaired Hessians, same masses (H 1.008) on both sides; modes matched by rank after sorting.")
    w("  Geometry shifts in Angstrom (CNEO - DFT).  'q' = the quantum proton(s) of the CNEO run (nuclear component index);")
    w("  where only one proton was quantum, the shift of the other H-X stretch is not a nuclear quantum effect.")
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
                qidx = [int(k[1:]) for k in c["qnuc"]]
                qlab = ", ".join((f"H{i} ({'donor' if i == dc.get('hd') else 'acceptor'})") for i in qidx) or "?"
                w("=" * 110)
                w(f"{start} {xc} {bas}   DFT bonding {dr.get('label', '?')}, CNEO bonding {dc.get('label', '?')}; q = {qlab}")
                w(f"  DFT  vibrations: {fmt(r['resR']['vib'])}")
                w(f"  CNEO vibrations: {fmt(c['resR']['vib'])}")
                if len(r["resR"]["vib"]) == len(c["resR"]["vib"]):
                    w("  CNEO - DFT     : " + ", ".join(f"{b - a:+.1f}" for a, b in zip(r["resR"]["vib"], c["resR"]["vib"])))
                if dr["valid"] and dc["valid"]:
                    w(f"  CNEO - DFT geometry: d(X-Hd) {dc['rXHd'] - dr['rXHd']:+.4f}  d(H..X) {dc['rHX'] - dr['rHX']:+.4f}  "
                      f"d(X-Ha) {dc['rXHa'] - dr['rXHa']:+.4f}  d(X..X) {dc['rXX'] - dr['rXX']:+.4f}  "
                      f"d(XH..X angle) {dc['a_XHX'] - dr['a_XHX']:+.1f}  d(tilt) {dc['a_HXH'] - dr['a_HXH']:+.1f} deg")
                w(f"  CNEO - DFT energy (Eh): {c['e_sp'] - r['e_sp']:+.6f}   (includes the quantum proton's zero-point-like energy; not a binding-energy shift)")

text = "\n".join(out) + "\n"
open(os.path.join(HERE, f"h{X.lower()}_h{X.lower()}_dft_classification.txt"), "w").write(text)
print(text)
