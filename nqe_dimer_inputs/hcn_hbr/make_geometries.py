"""Extract the final HCN...HBr geometries of the February 2026 CNEO-DFT runs (aug-cc-pV5Z excluded)
into the starting-geometry file read by make_inputs.py.

    python make_geometries.py <log dir> <out xyz>
"""
import glob
import math
import os
import re
import sys

LOGDIR = sys.argv[1]
XYZOUT = sys.argv[2]

BASIS = {"bas_two": "aug-cc-pVDZ", "bas_three": "aug-cc-pVTZ", "bas_four": "aug-cc-pVQZ"}
FUNC = {"func_one": "PBE", "func_two": "PW91", "func_three": "BP86", "func_four": "BLYP", "func_five": "B97"}
FUNC_ORDER = ["PBE", "PW91", "BP86", "BLYP", "B97"]
BAS_ORDER = ["aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]
STARTS = ["global_min", "local_min", "xlocal_min"]
ATOMS = ["H", "C", "N", "H", "Br"]
REAL = re.compile(r"[-+]?\d+\.\d+(?:e[-+]?\d+)?|[-+]?\d+e[-+]?\d+")
CYC = re.compile(r"cycle (\d+): E = (\S+)\s+dE = (\S+)\s+norm\(grad\) = (\S+)\n((?:[A-Z][a-z]?\s+\S+\s+\S+\s+\S+\n){5})")


def get_geom(text):
    for label, pat in (("summary", r"FINAL OPTIMIZED GEOMETRY: \n(\[\[.*?\]\])"),
                       ("post-opt array", r"BEGIN OPTIMIZED ATOMIC COORDINATES:\n-+\n(\[\[.*?\]\])")):
        m = re.search(pat, text, re.S)
        if m:
            nums = [float(x) for x in REAL.findall(m.group(1))]
            if len(nums) == 15:
                return [nums[i:i + 3] for i in range(0, 15, 3)], label
    cyc = CYC.findall(text)
    if cyc:
        geo = cyc[-1][4]
        xyz = [[float(x) for x in l.split()[1:]] for l in geo.strip().splitlines()]
        return xyz, "last optimizer step (5 decimals)"
    last = re.findall(r"^cycle (\d+): E = ", text, re.M)
    if last:
        n = last[-1]
        m = re.search(r"Geometry optimization cycle " + n +
                      r"\nCartesian coordinates \(Angstrom\)\n.*\n((?:\s+[A-Z][a-z]?\s+\S+\s+\S+\s+\S+.*\n){5})", text)
        if m:
            xyz = [[float(x) for x in l.split()[1:4]] for l in m.group(1).strip().splitlines()]
            return xyz, f"the last completed optimizer step {n} (6 decimals)"
    return None, None


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


records = {}
for path in sorted(glob.glob(os.path.join(LOGDIR, "*.log"))):
    parts = os.path.basename(path).split(".")
    bas, func, start, date = parts[1], parts[3], parts[4], parts[5]
    if bas not in BASIS or start not in STARTS:
        continue
    text = open(path, errors="replace").read()
    cycles = re.findall(r"cycle (\d+): E = (\S+)\s+dE = (\S+)\s+norm\(grad\) = (\S+)", text)
    if not cycles:
        print(f"skip (no optimizer cycles): {os.path.basename(path)}")
        continue
    xyz, src = get_geom(text)
    if xyz is None:
        print(f"skip (no geometry): {os.path.basename(path)}")
        continue
    opt_end = "END GEOMETRIC GEOMETRY OPTIMIZATION" in text
    capped = "Geometry optimization is not converged" in text
    ncyc = len(cycles)
    energy = float(cycles[-1][1])
    if opt_end and not capped:
        status = f"optimizer converged in {ncyc} steps"
    elif opt_end and capped:
        status = f"NOT converged: hit {ncyc}-step cap, Hessian run anyway"
    else:
        status = f"NOT finished: killed after {ncyc} steps"
    key = (start, BASIS[bas], FUNC[func])
    rec = dict(date=date, xyz=xyz, src=src, status=status, energy=energy, finished=opt_end, ncyc=ncyc,
               file=os.path.basename(path))
    old = records.get(key)
    if old is None or (rec["finished"], rec["ncyc"]) > (old["finished"], old["ncyc"]):
        if old is not None:
            print(f"duplicate for {key}: keeping {rec['file']} over {old['file']}")
        records[key] = rec
    elif old is not None:
        print(f"duplicate for {key}: keeping {old['file']} over {rec['file']}")

lines = []
missing = []
print(f"\n{'start':<11}{'xc':<6}{'basis':<13}{'E (Eh)':>16}  {'N..Br':>6} {'N..H(Br)':>8} {'H(C)..Br':>8}  status / source")
for start in STARTS:
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            rec = records.get((start, bas, xc))
            if rec is None:
                missing.append((start, bas, xc))
                continue
            xyz = rec["xyz"]
            n_br = dist(xyz[2], xyz[4])
            n_h = dist(xyz[2], xyz[3])
            hc_br = dist(xyz[0], xyz[4])
            note = rec["status"]
            if rec["src"] != "summary":
                note += f"; geometry from {rec['src']}"
            print(f"{start:<11}{xc:<6}{bas:<13}{rec['energy']:>16.9f}  {n_br:6.3f} {n_h:8.3f} {hc_br:8.3f}  {note}")
            lines.append("5")
            lines.append(f"{xc}/{bas} {start} E={rec['energy']:.9f} Eh {note}")
            for sym, (x, y, z) in zip(ATOMS, xyz):
                lines.append(f"{sym:<3}{x:16.8f}{y:16.8f}{z:16.8f}")

os.makedirs(os.path.dirname(XYZOUT), exist_ok=True)
open(XYZOUT, "w").write("\n".join(lines) + "\n")
print(f"\nwrote {len(lines) // 7} geometries to {XYZOUT}")
if missing:
    print("MISSING:", missing)
