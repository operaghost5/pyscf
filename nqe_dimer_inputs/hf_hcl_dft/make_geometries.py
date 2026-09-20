"""Starting geometries for the HF/HCl conventional-DFT re-optimization set (aug-cc-pV5Z excluded).

    python make_geometries.py <April log dir> <out xyz>

Both starts are the final geometries of the April 2026 conventional-DFT runs of the same functional
and basis, labelled by which molecule donates the hydrogen bond (read from the distances, not from
the April start label), the convention of ../hf_hcl/:

  hf_donor    F-H...Cl-H, HF the donor, HCl roughly perpendicular   (April 'global_min' runs)
  hcl_donor   Cl-H...F-H, HCl the donor, HF tilted                  (April 'local_min' runs)

Only the logs whose third file-name token is 'dft' are read.  Atom order everywhere: H0 F1 H2 Cl3
(H0-F1 = HF, H2-Cl3 = HCl); if a log has the hydrogens the other way round they are swapped and the
comment line says so.
"""
import glob
import math
import os
import re
import sys

LOGDIR = sys.argv[1]
XYZOUT = sys.argv[2]

METHOD_TOKEN = "dft"
BASIS = {"bas_two": "aug-cc-pVDZ", "bas_three": "aug-cc-pVTZ", "bas_four": "aug-cc-pVQZ"}
FUNC = {"func_one": "PBE", "func_two": "PW91", "func_three": "BP86", "func_four": "BLYP", "func_five": "B97"}
FUNC_ORDER = ["PBE", "PW91", "BP86", "BLYP", "B97"]
BAS_ORDER = ["aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]
ATOMS = ["H", "F", "H", "Cl"]
LABEL = {"F": "hf_donor", "Cl": "hcl_donor"}
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


def date_key(d):
    try:
        dd, mm, yy = d.split("-")
        return (int(yy), int(mm), int(dd))
    except ValueError:
        return (0, 0, 0)


def canonical(xyz):
    """Return (xyz with H0 on F1 and H2 on Cl3, donor symbol, swapped?)."""
    x = [row[:] for row in xyz]
    swapped = False
    if dist(x[0], x[1]) > dist(x[2], x[1]):          # H2 is the one on fluorine: swap the hydrogens
        x[0], x[2] = x[2], x[0]
        swapped = True
    hf_donates = dist(x[0], x[3]) < dist(x[2], x[1])   # H0...Cl shorter than H2...F
    return x, ("F" if hf_donates else "Cl"), swapped


records = {}
for path in sorted(glob.glob(os.path.join(LOGDIR, "*.log"))):
    parts = os.path.basename(path).split(".")
    if len(parts) < 6 or parts[2] != METHOD_TOKEN:
        continue
    bas, func, start, date = parts[1], parts[3], parts[4], parts[5]
    if bas not in BASIS or start not in ("global_min", "local_min"):
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
    if opt_end and not capped:
        status = f"optimizer converged in {len(cycles)} steps"
    elif opt_end:
        status = f"NOT converged: hit {len(cycles)}-step cap, Hessian run anyway"
    else:
        status = f"NOT finished: killed after {len(cycles)} steps"
    xyz, donor, swapped = canonical(xyz)
    if swapped:
        status += "; hydrogens swapped so that H0-F1 and H2-Cl3 are the covalent bonds"
    key = (LABEL[donor], BASIS[bas], FUNC[func])
    rec = dict(xyz=xyz, energy=energy, status=status, date=date, april=start, finished=opt_end,
               file=os.path.basename(path))
    old = records.get(key)
    if old is None or (rec["finished"], date_key(date)) > (old["finished"], date_key(old["date"])):
        if old is not None:
            print(f"two runs ended as {key}: keeping {rec['file']} over {old['file']}")
        records[key] = rec

lines = []
missing = []
print(f"{'start':<10}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'H0-F1':>6} {'H2-Cl3':>6} {'H0..Cl':>6} {'H2..F':>6} {'F..Cl':>6}  note")
for start in ("hf_donor", "hcl_donor"):
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            rec = records.get((start, bas, xc))
            if rec is None:
                missing.append((start, bas, xc))
                continue
            xyz = rec["xyz"]
            comment = (f"{xc}/{bas} {start} E={rec['energy']:.9f} Eh {rec['status']} "
                       f"(April 2026 conventional-DFT '{rec['april']}' run)")
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
if missing:
    print("MISSING (no April run ended as this isomer):", missing)
