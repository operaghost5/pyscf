"""Starting geometries for the (HBr)2 conventional-DFT re-optimization set (aug-cc-pV5Z excluded).

    python make_geometries.py <April log dir> <out xyz>

Two starts per functional and basis:

  global_min   the final geometry of the April 2026 conventional-DFT run of the same functional
               and basis (bent Cs hydrogen-bonded dimer; those runs started from an Angstrom
               geometry read as Bohr, but all 15 converged to the same structure and all are
               minima);
  c2h_saddle   a constructed C2h structure, the donor-acceptor interchange saddle point:
               Br...Br = 4.15 A, both H-Br bonds at 50 deg to the Br...Br axis in a trans
               arrangement, r(H-Br) taken from the acceptor H-Br of the April DFT geometry.
               The structure has exact inversion symmetry to floating-point precision,
               so an unconstrained minimizer stays on the C2h surface and converges to
               the saddle point; the input reports whether the symmetry survived.

Only the logs whose third file-name token is 'dft' are read.  Atom order everywhere:
H0 (donor H), Br1 (acceptor Br), H2 (acceptor H), Br3 (donor Br); H0-Br3 and Br1-H2 are the
covalent bonds, H0...Br1 the hydrogen bond.  If an April geometry has the donor at H2 the two
hydrogens are swapped and the comment line says so, the convention of ../hbr_hbr/.
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
ATOMS = ["H", "Br", "H", "Br"]
NATM = 4
REAL = re.compile(r"[-+]?\d+\.\d+(?:e[-+]?\d+)?|[-+]?\d+e[-+]?\d+")

C2H_FF = 4.15          # Angstrom, Br...Br distance of the constructed saddle start
C2H_TILT_DEG = 50.0    # angle between each H-Br bond and the Br...Br axis


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


def c2h_start(r_hf, R=C2H_FF, tilt_deg=C2H_TILT_DEG):
    """H0 Br1 H2 Br3 coordinates (Angstrom) of the trans-tilted C2h structure, inversion-symmetric."""
    th = math.radians(tilt_deg)
    fa = [-0.5 * R, 0.0, 0.0]                                   # donor Br (atom 3)
    ha = [fa[0] + r_hf * math.cos(th), r_hf * math.sin(th), 0.0]  # donor H (atom 0), pointing at Br1
    fb = [-fa[0], -fa[1], -fa[2]]                                # acceptor Br (atom 1) = inversion of Br3
    hb = [-ha[0], -ha[1], -ha[2]]                                # acceptor H (atom 2) = inversion of H0
    return [ha, fb, hb, fa]


def donor_first(xyz):
    """Return (xyz with the donor H at index 0, swapped?).  H0-Br3 / H2-Br1 covalent, H0...Br1 the H bond."""
    d = dist
    # covalent partner of each H: the nearer Br
    p0 = 3 if d(xyz[0], xyz[3]) < d(xyz[0], xyz[1]) else 1
    p2 = 3 if d(xyz[2], xyz[3]) < d(xyz[2], xyz[1]) else 1
    if p0 == p2:
        raise ValueError("both hydrogens bonded to the same bromine")
    other = {1: 3, 3: 1}
    hb0, hb2 = d(xyz[0], xyz[other[p0]]), d(xyz[2], xyz[other[p2]])
    donor = 0 if hb0 < hb2 else 2
    x = [row[:] for row in xyz]
    swapped = False
    if donor == 2:                                   # put the donor H at index 0
        x[0], x[2] = x[2], x[0]
        swapped = True
        p0, p2 = p2, p0
    if p0 == 1:                                      # donor H bonded to Br1: swap the bromines so that Br3 is the donor Br
        x[1], x[3] = x[3], x[1]
        swapped = True
    return x, swapped


records = {}
for path in sorted(glob.glob(os.path.join(LOGDIR, "*.log"))):
    parts = os.path.basename(path).split(".")
    if len(parts) < 6 or parts[2] != METHOD_TOKEN:
        continue
    bas, func, start, date = parts[1], parts[3], parts[4], parts[5]
    if bas not in BASIS or start != "global_min":
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
    xyz, swapped = donor_first(xyz)
    if swapped:
        status += "; hydrogens/bromines reordered so that H0-Br3 donates to Br1"
    key = (BASIS[bas], FUNC[func])
    rec = dict(xyz=xyz, energy=energy, status=status, date=date, finished=opt_end, file=os.path.basename(path))
    old = records.get(key)
    if old is None or (rec["finished"], date_key(date)) > (old["finished"], date_key(old["date"])):
        if old is not None:
            print(f"duplicate for {key}: keeping {rec['file']} over {old['file']}")
        records[key] = rec

lines = []
missing = []
print(f"{'start':<11}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'Br3-H0':>6} {'H0..Br1':>7} {'Br1-H2':>6} {'Br1..Br3':>7} {'H2..Br3':>7}  note")
for start in ("global_min", "c2h_saddle"):
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            rec = records.get((bas, xc))
            if rec is None:
                missing.append((start, bas, xc))
                continue
            if start == "global_min":
                xyz = rec["xyz"]
                comment = (f"{xc}/{bas} global_min E={rec['energy']:.9f} Eh {rec['status']} "
                           f"(April 2026 conventional-DFT run, Angstrom start read as Bohr)")
            else:
                r_hf = dist(rec["xyz"][1], rec["xyz"][2])            # April DFT acceptor H-Br length
                xyz = c2h_start(r_hf)
                comment = (f"{xc}/{bas} c2h_saddle E=n/a Eh constructed C2h interchange-saddle start: "
                           f"Br..Br {C2H_FF:.2f} A, r(H-Br) {r_hf:.4f} A from the April DFT acceptor bond, "
                           f"tilt {C2H_TILT_DEG:.0f} deg, exact inversion symmetry")
            d = dist
            print(f"{start:<11}{xc:<6}{bas:<13}{rec['energy'] if start == 'global_min' else float('nan'):>16.9f} "
                  f"{d(xyz[3], xyz[0]):6.3f} {d(xyz[0], xyz[1]):7.3f} {d(xyz[1], xyz[2]):6.3f} "
                  f"{d(xyz[1], xyz[3]):7.3f} {d(xyz[2], xyz[3]):7.3f}  {comment.split(' Eh ', 1)[1][:70]}")
            lines.append(str(NATM))
            lines.append(comment)
            for sym, (x, y, z) in zip(ATOMS, xyz):
                lines.append(f"{sym:<3}{x:16.8f}{y:16.8f}{z:16.8f}")

os.makedirs(os.path.dirname(XYZOUT), exist_ok=True)
open(XYZOUT, "w").write("\n".join(lines) + "\n")
print(f"\nwrote {len(lines) // (NATM + 2)} geometries to {XYZOUT}")
if missing:
    print("MISSING:", missing)
