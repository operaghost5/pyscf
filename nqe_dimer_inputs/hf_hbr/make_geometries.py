"""Starting geometries for the HF/HBr re-optimization set (aug-cc-pV5Z excluded).

    python make_geometries.py <April log dir> <out xyz>

Both starts are the final geometries of the April 2026 CNEO-DFT runs of the same functional and
basis, relabelled by which molecule donates the hydrogen bond, because the April energies cannot
rank the two isomers (each run had only its own donor proton quantum):

  hf_donor    April 'global_min':  F-H...Br-H, HF the donor, HBr roughly perpendicular
  hbr_donor   April 'local_min':   Br-H...F-H, HBr the donor, HF tilted

Atom order everywhere: H0 F1 H2 Br3 (H0-F1 = HF, H2-Br3 = HBr).
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
RELABEL = {"global_min": "hf_donor", "local_min": "hbr_donor"}
ATOMS = ["H", "F", "H", "Br"]
NATM = 4
REAL = re.compile(r"[-+]?\d+\.\d+(?:e[-+]?\d+)?|[-+]?\d+e[-+]?\d+")


def get_geom(text):
    for label, pat in (("summary", r"FINAL OPTIMIZED GEOMETRY: \n(\[\[.*?\]\])"),
                       ("post-opt array", r"BEGIN OPTIMIZED ATOMIC COORDINATES:\n-+\n(\[\[.*?\]\])")):
        m = re.search(pat, text, re.S)
        if m:
            nums = [float(x) for x in REAL.findall(m.group(1))]
            if len(nums) == 3 * NATM:
                return [nums[i:i + 3] for i in range(0, 3 * NATM, 3)], label
    return None, None


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


records = {}
for path in sorted(glob.glob(os.path.join(LOGDIR, "*.log"))):
    parts = os.path.basename(path).split(".")
    bas, func, start, date = parts[1], parts[3], parts[4], parts[5]
    if bas not in BASIS or start not in RELABEL:
        continue
    text = open(path, errors="replace").read()
    cycles = re.findall(r"^cycle (\d+): E = (\S+)\s+dE = (\S+)\s+norm\(grad\) = (\S+)", text, re.M)
    xyz, src = get_geom(text)
    if not cycles or xyz is None:
        print(f"skip: {os.path.basename(path)}")
        continue
    opt_end = "END GEOMETRIC GEOMETRY OPTIMIZATION" in text
    capped = "Geometry optimization is not converged" in text
    m = re.search(r"FINAL SINGLE POINT ENERGY:\s+(\S+)", text)
    energy = float(m.group(1)) if m else float(cycles[-1][1])
    qn = ",".join(sorted(set(re.findall(r"\|g_(n\d+)\|=", text))))
    if opt_end and not capped:
        status = f"optimizer converged in {len(cycles)} steps"
    elif opt_end:
        status = f"NOT converged: hit {len(cycles)}-step cap, Hessian run anyway"
    else:
        status = f"NOT finished: killed after {len(cycles)} steps"
    records[(RELABEL[start], BASIS[bas], FUNC[func])] = dict(xyz=xyz, energy=energy, status=status, date=date,
                                                              april=start, qn=qn)

lines = []
print(f"{'start':<10}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'H0-F1':>6} {'H2-Br3':>6} {'H0..Br':>6} {'H2..F':>6} {'F..Br':>6}  note")
for start in ("hf_donor", "hbr_donor"):
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            rec = records[(start, bas, xc)]
            xyz = rec["xyz"]
            comment = (f"{xc}/{bas} {start} E={rec['energy']:.9f} Eh {rec['status']} "
                       f"(April 2026 '{rec['april']}' run, quantum nucleus {rec['qn']} only, i.e. the donor proton)")
            d = dist
            print(f"{start:<10}{xc:<6}{bas:<13}{rec['energy']:>16.9f} {d(xyz[0], xyz[1]):6.3f} {d(xyz[2], xyz[3]):6.3f} "
                  f"{d(xyz[0], xyz[3]):6.3f} {d(xyz[2], xyz[1]):6.3f} {d(xyz[1], xyz[3]):6.3f}  {rec['status']}")
            lines.append(str(NATM))
            lines.append(comment)
            for sym, (x, y, z) in zip(ATOMS, xyz):
                lines.append(f"{sym:<3}{x:16.8f}{y:16.8f}{z:16.8f}")

os.makedirs(os.path.dirname(XYZOUT), exist_ok=True)
open(XYZOUT, "w").write("\n".join(lines) + "\n")
print(f"\nwrote {len(lines) // (NATM + 2)} geometries to {XYZOUT}")
