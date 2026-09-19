"""HCl/HBr dimer: classification of the April 2026 CNEO-DFT runs (aug-cc-pV5Z excluded).

Reuses the parsing / projection functions of ../hcncl/reanalyze_hessians2.py and adds
the translational sum rule per atom, the diagonal-block repaired spectrum, the geometry
descriptors and the optimizer status.  Pure Python.
"""
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, "reanalyze_hessians2.py")).read()   # generic parser/projector, copied alongside
exec(SRC.split("\nlogs = sorted(")[0])          # functions and constants only, no main loop

LOGDIR = sys.argv[1]
BASNAME = {"bas_two": "aug-cc-pVDZ", "bas_three": "aug-cc-pVTZ", "bas_four": "aug-cc-pVQZ"}
FUNC_ORDER = ["PBE", "PW91", "BP86", "BLYP", "B97"]
BAS_ORDER = ["aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]


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


records = {}
for p in sorted(glob.glob(os.path.join(LOGDIR, "*.log"))):
    parts = os.path.basename(p).split(".")
    bas, func, start, date = parts[1], parts[3], parts[4], parts[5]
    if bas not in BASNAME:
        continue
    t = open(p, errors="replace").read()
    cycles = re.findall(r"^cycle (\d+): E = (\S+)\s+dE = (\S+)\s+norm\(grad\) = (\S+)", t, re.M)
    rec = dict(start=start, xc=FUNC[func], basis=BASNAME[bas], date=date, ncyc=len(cycles),
               lines=t.count("\n"), scf_started="BEGIN INITIAL SCF RUN" in t,
               opt_end="END GEOMETRIC GEOMETRY OPTIMIZATION" in t,
               capped="Geometry optimization is not converged" in t,
               hess_end="END HESSIAN CALCULATION" in t, last_line=t.strip().splitlines()[-1] if t.strip() else "")
    m = re.search(r"FINAL SINGLE POINT ENERGY:\s+(\S+)", t)
    rec["e_sp"] = float(m.group(1)) if m else None
    rec["qnuc"] = sorted(set(re.findall(r"\|g_(n\d+)\|=", t)))
    if rec["hess_end"]:
        symbols = parse_symbols(t)
        natm = len(symbols)
        H = parse_hessian(t, natm)
        g = parse_geom(t, natm)
        res = analyse(H, g, symbols)
        viol = sum_rule(H, natm)
        resR = analyse(repair(H, natm), g, symbols)
        logged = re.search(r"BEGIN HARMONIC ANALYSIS.*?\n-+\n(.*?)\n-+\nEND HARMONIC", t, re.S)
        rec.update(symbols=symbols, xyz=g, res=res, resR=resR, viol=viol,
                   logged=[float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:e[-+]?\d+)?", logged.group(1))] if logged else None)
    records[(start, rec["basis"], rec["xc"])] = rec

out = []
w = out.append
w("HCl/HBr dimer, CNEO-DFT April 2026 runs: classification (aug-cc-pV5Z excluded)")
w("Masses: PySCF isotope-averaged (H 1.008, Cl 35.45, Br 79.904). Atom order in every log: H0 Cl1 H2 Br3; which H belongs to which molecule is read from the distances.")
w("")
w("PART 0: run status")
w(f"{'start':<11}{'xc':<6}{'basis':<13}{'opt steps':>9}  {'quantum nuclei':<16}status")
for start in ("global_min", "local_min"):
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None:
                w(f"{start:<11}{xc:<6}{bas:<13}{'-':>9}  no log")
                continue
            if not r["scf_started"]:
                st = f"DIED BEFORE THE FIRST SCF ({r['lines']} lines; last line: '{r['last_line']}')"
            elif not r["opt_end"]:
                st = "optimization did not finish"
            elif r["capped"]:
                st = "optimizer hit its step cap; Hessian run anyway" if r["hess_end"] else "optimizer hit its step cap; no Hessian"
            else:
                st = "optimizer converged; Hessian done" if r["hess_end"] else "optimizer converged; no Hessian"
            w(f"{start:<11}{xc:<6}{bas:<13}{r['ncyc']:>9d}  {','.join(r['qnuc']) or '?':<16}{st}")
w("")
w("PART 1: geometry, energy and classification of the runs that produced a Hessian")
w("  bonding read from distances: hCl = the H bonded to Cl, hBr = the H bonded to Br. r(hCl..Br) = Cl-H...Br hydrogen bond (HCl donor),")
w("  r(hBr..Cl) = Br-H...Cl (HBr donor); donor = the shorter of the two; angle = X-H...Y of that bond; tilt = H...Y-H at the acceptor;")
w("  'q' = which molecule's proton was quantum in the run (from the nuclear component index).")
w(f"{'start':<11}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'H-Cl':>6} {'H-Br':>6} {'hCl.Br':>6} {'hBr.Cl':>6} {'Cl..Br':>6} {'donor':>6} {'XH..Y':>6} {'tilt':>6} {'q':>4}  "
  f"{'Cl viol':>8} {'Br viol':>8}  vibrations (cm^-1), 3N-6 = 6            imag  class")
class_summary = []
for start in ("global_min", "local_min"):
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None or not r.get("hess_end"):
                continue
            x = r["xyz"]
            F, Cl = 1, 3                                              # here: X = Cl (index 1), Y = Br (index 3)
            hF = 0 if dist(x[0], x[F]) < dist(x[2], x[F]) else 2      # the H covalently bound to Cl
            hCl = 2 if hF == 0 else 0                                 # the H covalently bound to Br
            rHF, rHCl = dist(x[hF], x[F]), dist(x[hCl], x[Cl])
            rHFCl, rHClF, rFCl = dist(x[hF], x[Cl]), dist(x[hCl], x[F]), dist(x[F], x[Cl])
            if rHFCl < rHClF:
                donor, a1, a2 = "HCl", angle(x[F], x[hF], x[Cl]), angle(x[hF], x[Cl], x[hCl])
            else:
                donor, a1, a2 = "HBr", angle(x[Cl], x[hCl], x[F]), angle(x[hCl], x[F], x[hF])
            qidx = [int(k[1:]) for k in r["qnuc"]]
            q = "+".join(("HCl" if i == hF else "HBr") for i in qidx) or "?"
            vib = r["res"]["vib"]
            nimag = sum(1 for f in vib if f < 0)
            cls = "MINIMUM" if nimag == 0 else f"SADDLE (order {nimag})"
            w(f"{start:<11}{xc:<6}{bas:<13}{r['e_sp']:>16.9f} {rHF:6.3f} {rHCl:6.3f} {rHFCl:6.3f} {rHClF:6.3f} {rFCl:6.3f} {donor:>6} {a1:6.1f} {a2:6.1f} {q:>4}  "
              f"{r['viol'][1]:8.1e} {r['viol'][3]:8.1e}  {fmt(vib):<40s} {nimag:>3d}  {cls}")
            class_summary.append((start, bas, xc, cls, r['e_sp'], q))
w("")
w("PART 2: external modes and the halogen sum-rule defect")
w("  'full 3N' is the unprojected spectrum: 6 external modes should be ~0 but are displaced by the")
w("  bromine and chlorine diagonal-block defects; 'repaired' resets each diagonal block to minus the sum of its")
w("  off-diagonal blocks (the production repair).  The vibrations are compared before/after.")
for start in ("global_min", "local_min"):
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            r = records.get((start, bas, xc))
            if r is None or not r.get("hess_end"):
                continue
            v = r["viol"]
            w("=" * 100)
            w(f"{start} {xc} {bas}   sum rule per atom (Eh/Bohr^2): "
              + ", ".join(f"{s}{i} {vv:.1e}" for i, (s, vv) in enumerate(zip(r['symbols'], v))))
            w("  PySCF logged (3N-6)      : " + fmt(r["logged"]))
            w("  re-diagonalized, projected: " + fmt(r["res"]["vib"]))
            w("  full 3N, raw              : " + fmt(r["res"]["full"]))
            w("  full 3N, repaired         : " + fmt(r["resR"]["full"]))
            w("  repaired, projected       : " + fmt(r["resR"]["vib"]))
            dmax = max(abs(a - b) for a, b in zip(r["res"]["vib"], r["resR"]["vib"]))
            w(f"  max |raw - repaired| over the 6 vibrations: {dmax:.2f} cm^-1")
w("")
w("PART 3: classification on the RAW and on the REPAIRED Hessian (the repaired one is the physical answer:")
w("  a halogen diagonal-block defect too large for the TR projection to absorb leaks a displaced translation into")
w("  the 'vibrations'; the repaired spectra are the ones to compare across bases)")
w(f"{'start':<11}{'xc':<6}{'basis':<13} {'raw class':<18}{'repaired class':<16}  repaired vibrations (cm^-1)")
for start in ("global_min", "local_min"):
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
K = 627.5095
w("  energy of the local_min run relative to the global_min run of the same functional/basis (kcal/mol), with the")
w("  quantum proton of each run (see column q above): comparable only if the same molecule's proton was quantum in both.")
for bas in BAS_ORDER:
    for xc in FUNC_ORDER:
        eg = [c[4] for c in class_summary if c[0] == "global_min" and c[1] == bas and c[2] == xc]
        el = [c[4] for c in class_summary if c[0] == "local_min" and c[1] == bas and c[2] == xc]
        qg = [c[5] for c in class_summary if c[0] == "global_min" and c[1] == bas and c[2] == xc]
        ql = [c[5] for c in class_summary if c[0] == "local_min" and c[1] == bas and c[2] == xc]
        if eg and el:
            w(f"    {xc:<6}{bas:<13}{(el[0] - eg[0]) * K:+8.2f}   quantum proton: global {qg[0]}, local {ql[0]}")
text = "\n".join(out) + "\n"
open(os.path.join(HERE, "hcl_hbr_classification.txt"), "w").write(text)
print(text)
