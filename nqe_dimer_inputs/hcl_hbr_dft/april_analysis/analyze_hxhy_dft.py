"""HX/HY heterodimers (HF/HBr, HF/HCl, HCl/HBr): classification of the April 2026 conventional-DFT runs
(aug-cc-pV5Z excluded).

Reuses the parsing / projection functions of reanalyze_hessians2.py (copied alongside) and
adds the translational sum rule per atom, the diagonal-block repaired spectrum, the geometry
descriptors (covalent partners and the donor read from the distances, so the isomer is named
by its donor whatever the run was called), the optimizer status, the isomer energy ordering
(directly comparable here because no nucleus is quantum) and, when a directory of CNEO-DFT logs
is given as the second argument, the CNEO - DFT frequency and geometry shifts for the runs that
have a Hessian in both, with the quantum proton(s) of each CNEO run noted.

    python analyze_hxhy_dft.py <DFT log dir> [<CNEO-DFT log dir>]

Atom order in every log: H0 X1 H2 Y3 with X the lighter halogen; the halogens are read from the
logs.  Pure Python.
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


def descriptors(x, symbols):
    """H0 X1 H2 Y3 (any order of the two halogens).  Covalent partners and the donor from the distances."""
    hal = [i for i, s in enumerate(symbols) if s != "H"]
    hyd = [i for i, s in enumerate(symbols) if s == "H"]
    partner = {h: min(hal, key=lambda a: dist(x[h], x[a])) for h in hyd}
    if len(set(partner.values())) < 2:
        return dict(valid=False)
    other = {hal[0]: hal[1], hal[1]: hal[0]}
    hb = {h: dist(x[h], x[other[partner[h]]]) for h in hyd}
    hd = min(hb, key=hb.get)
    ha = [h for h in hyd if h != hd][0]
    Xd, Xa = partner[hd], partner[ha]
    d = dict(valid=True, hd=hd, ha=ha, Xd=Xd, Xa=Xa, donor=f"H{symbols[Xd]}",
             rXHd=dist(x[hd], x[Xd]), rHX=hb[hd], rXHa=dist(x[ha], x[Xa]), rXY=dist(x[Xd], x[Xa]),
             rHXb=hb[ha], a_XHY=angle(x[Xd], x[hd], x[Xa]), a_HYH=angle(x[hd], x[Xa], x[ha]))
    d["label"] = f"{symbols[Xd]}{Xd}-H{hd}...{symbols[Xa]}{Xa}-H{ha}"
    d["bound"] = d["rHX"] < 4.0
    # covalent bonds per molecule, keyed by halogen symbol, for the CNEO - DFT geometry comparison
    d["bond"] = {symbols[partner[h]]: dist(x[h], x[partner[h]]) for h in hyd}
    return d


records, dups = read_logs(LOGDIR, "dft")
cneo, cdups = read_logs(CNEODIR, "cneodft") if CNEODIR else ({}, [])
symbols0 = None
for r in records.values():
    if r.get("symbols"):
        symbols0 = r["symbols"]
        break
if symbols0 is None:
    sys.exit("no DFT log with a Hessian found")
HAL = [s for s in symbols0 if s != "H"]
X, Y = HAL[0], HAL[1]
NAME = f"H{X}/H{Y}"

out = []
w = out.append
w(f"{NAME}, conventional DFT (all nuclei classical), April 2026 runs: classification (aug-cc-pV5Z excluded)")
w(f"Masses: PySCF isotope-averaged (H 1.008, {X} {MASS[X]}, {Y} {MASS[Y]}).  Atom order in every log: {' '.join(f'{s}{i}' for i, s in enumerate(symbols0))};")
w("the covalent partners and the donor are read from the distances, so each run is named by its donor whatever its start label.")
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
w("  bonding = donor halogen and H, acceptor halogen and H, from the distances; X-Hd = donor covalent bond, H..Y = hydrogen bond,")
w("  Y-Ha = acceptor covalent bond, X..Y; angles X-H...Y (linearity of the H bond) and H...Y-H (tilt of the acceptor);")
w(f"  '{X} viol' / '{Y} viol' = sum-rule violation of each halogen (Eh/Bohr^2).")
w(f"{'start':<11}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'bonding':<16} {'donor':>5} {'X-Hd':>6} {'H..Y':>6} {'Y-Ha':>6} {'X..Y':>6} {'XH..Y':>6} {'H..YH':>6} "
  f"{X+' viol':>8} {Y+' viol':>8}  vibrations (cm^-1), 3N-6 = 6                                      imag  class")
class_summary = []
for start in STARTS:
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None or not r.get("hess_end"):
                continue
            d = descriptors(r["xyz"], r["symbols"])
            vib = r["res"]["vib"]
            nimag = sum(1 for f in vib if f < 0)
            cls = "MINIMUM" if nimag == 0 else f"SADDLE (order {nimag})"
            hal_idx = [i for i, s in enumerate(r["symbols"]) if s != "H"]
            vX, vY = r["viol"][hal_idx[0]], r["viol"][hal_idx[1]]
            if not d["valid"]:
                w(f"{start:<11}{xc:<6}{bas:<13}{r['e_sp']:>16.9f} {'not an HX/HY dimer':<16}")
                continue
            flag = "" if d["bound"] else "  UNBOUND (H..Y > 4 A)"
            w(f"{start:<11}{xc:<6}{bas:<13}{r['e_sp']:>16.9f} {d['label']:<16} {d['donor']:>5} {d['rXHd']:6.3f} {d['rHX']:6.3f} {d['rXHa']:6.3f} {d['rXY']:6.3f} "
              f"{d['a_XHY']:6.1f} {d['a_HYH']:6.1f} {vX:8.1e} {vY:8.1e}  {fmt(vib):<58s} {nimag:>3d}  {cls}{flag}")
            class_summary.append((start, bas, xc, cls, r["e_sp"], d["donor"], d["bound"]))
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
w(f"  isomer energies: E({NAME.split('/')[1]}-donor run) - E({NAME.split('/')[0]}-donor run) in kcal/mol, same functional and basis; all nuclei are")
w("  classical, so the two isomers are directly comparable (the April CNEO energies were not: a different proton was quantum in each):")
w(f"    {'xc':<6}{'basis':<13}{'dE':>8}   donors (global_min run / local_min run)")
for bas in BAS_ORDER:
    for xc in FUNC_ORDER:
        rows = {c[0]: c for c in class_summary if c[1] == bas and c[2] == xc}
        if len(rows) < 2:
            continue
        by_donor = {c[5]: c for c in rows.values()}
        dX, dY = f"H{X}", f"H{Y}"
        if dX in by_donor and dY in by_donor:
            de = (by_donor[dY][4] - by_donor[dX][4]) * K
            w(f"    {xc:<6}{bas:<13}{de:>+8.2f}   {rows['global_min'][5]} / {rows['local_min'][5]}"
              + ("" if all(c[6] for c in rows.values()) else "   (an unbound structure is involved)"))
        else:
            w(f"    {xc:<6}{bas:<13}{'-':>8}   both runs ended as the {list(by_donor)[0]}-donor isomer")

if cneo:
    w("")
    w("PART 4: CNEO-DFT minus DFT, same functional, basis and start, for the runs with a Hessian in both sets")
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
                dr, dc = descriptors(r["xyz"], r["symbols"]), descriptors(c["xyz"], c["symbols"])
                qidx = [int(k[1:]) for k in c["qnuc"]]
                qlab = ", ".join(f"H{i} (H{c['symbols'][min([j for j, s in enumerate(c['symbols']) if s != 'H'], key=lambda a: dist(c['xyz'][i], c['xyz'][a]))]}"
                                 f", {'donor' if i == dc.get('hd') else 'acceptor'})" for i in qidx) or "?"
                w("=" * 110)
                w(f"{start} {xc} {bas}   DFT bonding {dr.get('label', '?')}, CNEO bonding {dc.get('label', '?')}; q = {qlab}")
                w(f"  DFT  vibrations: {fmt(r['resR']['vib'])}")
                w(f"  CNEO vibrations: {fmt(c['resR']['vib'])}")
                if len(r["resR"]["vib"]) == len(c["resR"]["vib"]):
                    w("  CNEO - DFT     : " + ", ".join(f"{b - a:+.1f}" for a, b in zip(r["resR"]["vib"], c["resR"]["vib"])))
                if dr["valid"] and dc["valid"] and dr["donor"] == dc["donor"]:
                    w(f"  CNEO - DFT geometry: d(H-{X}) {dc['bond'][X] - dr['bond'][X]:+.4f}  d(H-{Y}) {dc['bond'][Y] - dr['bond'][Y]:+.4f}  "
                      f"d(H..Y) {dc['rHX'] - dr['rHX']:+.4f}  d(X..Y) {dc['rXY'] - dr['rXY']:+.4f}  "
                      f"d(XH..Y angle) {dc['a_XHY'] - dr['a_XHY']:+.1f}  d(tilt) {dc['a_HYH'] - dr['a_HYH']:+.1f} deg")
                elif dr["valid"] and dc["valid"]:
                    w("  (different donors in the DFT and CNEO runs: no geometry comparison)")
                w(f"  CNEO - DFT energy (Eh): {c['e_sp'] - r['e_sp']:+.6f}   (includes the quantum proton's zero-point-like energy; not a binding-energy shift)")

text = "\n".join(out) + "\n"
open(os.path.join(HERE, f"h{X.lower()}_h{Y.lower()}_dft_classification.txt"), "w").write(text)
print(text)
