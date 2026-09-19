"""Starting geometries for the HCl/HBr re-optimization set (aug-cc-pV5Z excluded).

    python make_geometries.py <April log dir> <out xyz>

Both starts are the final geometries of the April 2026 CNEO-DFT runs of the same functional and
basis, relabelled by which molecule donates the hydrogen bond, because the April energies cannot
rank the two isomers (each run had only its own donor proton quantum):

  hcl_donor   April 'global_min':  Cl-H...Br-H, HCl the donor, HBr roughly perpendicular
  hbr_donor   April 'local_min':   Br-H...Cl-H, HBr the donor, HCl roughly perpendicular

Atom order in the output: H0 Cl1 H2 Br3 with H0-Cl1 = HCl and H2-Br3 = HBr in EVERY block (the two
hydrogens are swapped if a log has them the other way round; which H belongs to which molecule is
decided from the distances, never assumed). For a run that died before the final geometry was
printed (B97/aug-cc-pVQZ global_min, killed at step 438), the geometry of the last completed
optimizer step is used and the comment line says so.
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
RELABEL = {"global_min": "hcl_donor", "local_min": "hbr_donor"}
ATOMS = ["H", "Cl", "H", "Br"]
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
    last = re.findall(r"^cycle (\d+): E = ", text, re.M)
    if last:
        n = last[-1]
        m = re.search(r"Geometry optimization cycle " + n +
                      r"\nCartesian coordinates \(Angstrom\)\n.*\n((?:\s+[A-Z][a-z]?\s+\S+\s+\S+\s+\S+.*\n){" + str(NATM) + "})", text)
        if m:
            xyz = [[float(x) for x in l.split()[1:4]] for l in m.group(1).strip().splitlines()]
            return xyz, f"the last completed optimizer step {n}"
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
    swapped = False
    if dist(xyz[2], xyz[1]) < dist(xyz[0], xyz[1]):      # H2 is the one bonded to Cl: swap the hydrogens
        xyz = [xyz[2], xyz[1], xyz[0], xyz[3]]
        swapped = True
    assert dist(xyz[0], xyz[1]) < 1.6 and dist(xyz[2], xyz[3]) < 1.8, (path, xyz)
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
    key = (RELABEL[start], BASIS[bas], FUNC[func])
    rec = dict(xyz=xyz, energy=energy, status=status, date=date, april=start, qn=qn, swapped=swapped,
               src=src, finished=opt_end, ncyc=len(cycles))
    old = records.get(key)
    if old is None or (rec["finished"], rec["ncyc"]) > (old["finished"], old["ncyc"]):
        records[key] = rec

lines = []
print(f"{'start':<10}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'H0-Cl1':>6} {'H2-Br3':>6} {'H0..Br':>6} {'H2..Cl':>6} {'Cl..Br':>6}  note")
for start in ("hcl_donor", "hbr_donor"):
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            rec = records[(start, bas, xc)]
            xyz = rec["xyz"]
            qwhich = "the donor proton"                    # n0 = HCl proton in global_min, n2 = HBr proton in local_min
            comment = (f"{xc}/{bas} {start} E={rec['energy']:.9f} Eh {rec['status']} "
                       f"(April 2026 '{rec['april']}' run, quantum nucleus {rec['qn']} only = {qwhich}"
                       + ("; hydrogens swapped into the H0-Cl1 / H2-Br3 order" if rec['swapped'] else "")
                       + (f"; geometry from {rec['src']}" if rec['src'] != "summary" else "") + ")")
            d = dist
            print(f"{start:<10}{xc:<6}{bas:<13}{rec['energy']:>16.9f} {d(xyz[0], xyz[1]):6.3f} {d(xyz[2], xyz[3]):6.3f} "
                  f"{d(xyz[0], xyz[3]):6.3f} {d(xyz[2], xyz[1]):6.3f} {d(xyz[1], xyz[3]):6.3f}  {rec['status']}{' (H swapped)' if rec['swapped'] else ''}{'' if rec['src'] == 'summary' else ' [' + rec['src'] + ']'}")
            lines.append(str(NATM))
            lines.append(comment)
            for sym, (x, y, z) in zip(ATOMS, xyz):
                lines.append(f"{sym:<3}{x:16.8f}{y:16.8f}{z:16.8f}")

os.makedirs(os.path.dirname(XYZOUT), exist_ok=True)
open(XYZOUT, "w").write("\n".join(lines) + "\n")
print(f"\nwrote {len(lines) // (NATM + 2)} geometries to {XYZOUT}")
