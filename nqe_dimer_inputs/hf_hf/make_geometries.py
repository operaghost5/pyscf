"""Starting geometries for the (HF)2 re-optimization set (aug-cc-pV5Z excluded).

    python make_geometries.py <April log dir> <out xyz>

Two starts per functional and basis:

  global_min   the final geometry of the April 2026 CNEO-DFT run of the same functional
               and basis (bent Cs hydrogen-bonded dimer; those runs had only the donor
               proton quantum and started from an Angstrom geometry read as Bohr, but
               all 15 converged to the same structure and all are minima);
  c2h_saddle   a constructed C2h structure, the donor-acceptor interchange saddle point:
               F...F = 2.78 A, both H-F bonds at 55 deg to the F...F axis in a trans
               arrangement, r(H-F) taken from the acceptor H-F of the April geometry.
               The structure has exact inversion symmetry to floating-point precision,
               so an unconstrained minimizer stays on the C2h surface and converges to
               the saddle point; the input reports whether the symmetry survived.

Atom order everywhere: H0 (donor H), F1 (acceptor F), H2 (acceptor H), F3 (donor F);
H0-F3 and F1-H2 are the covalent bonds, H0...F1 the hydrogen bond.
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
ATOMS = ["H", "F", "H", "F"]
NATM = 4
REAL = re.compile(r"[-+]?\d+\.\d+(?:e[-+]?\d+)?|[-+]?\d+e[-+]?\d+")

C2H_FF = 2.78          # Angstrom, F...F distance of the constructed saddle start
C2H_TILT_DEG = 55.0    # angle between each H-F bond and the F...F axis


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


def c2h_start(r_hf, R=C2H_FF, tilt_deg=C2H_TILT_DEG):
    """H0 F1 H2 F3 coordinates (Angstrom) of the trans-tilted C2h structure, inversion-symmetric."""
    th = math.radians(tilt_deg)
    fa = [-0.5 * R, 0.0, 0.0]                                   # donor F (atom 3)
    ha = [fa[0] + r_hf * math.cos(th), r_hf * math.sin(th), 0.0]  # donor H (atom 0), pointing at F1
    fb = [-fa[0], -fa[1], -fa[2]]                                # acceptor F (atom 1) = inversion of F3
    hb = [-ha[0], -ha[1], -ha[2]]                                # acceptor H (atom 2) = inversion of H0
    return [ha, fb, hb, fa]


records = {}
for path in sorted(glob.glob(os.path.join(LOGDIR, "*.log"))):
    parts = os.path.basename(path).split(".")
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
    records[(BASIS[bas], FUNC[func])] = dict(xyz=xyz, energy=energy, status=status, date=date)

lines = []
print(f"{'start':<11}{'xc':<6}{'basis':<13}{'E (Eh)':>16} {'F3-H0':>6} {'H0..F1':>7} {'F1-H2':>6} {'F1..F3':>7} {'H2..F3':>7}  note")
for start in ("global_min", "c2h_saddle"):
    for bas in BAS_ORDER:
        for xc in FUNC_ORDER:
            rec = records[(bas, xc)]
            if start == "global_min":
                xyz = rec["xyz"]
                comment = (f"{xc}/{bas} global_min E={rec['energy']:.9f} Eh {rec['status']} "
                           f"(April 2026 run, donor proton quantum only, Angstrom start read as Bohr)")
            else:
                r_hf = dist(rec["xyz"][1], rec["xyz"][2])            # April acceptor H-F length
                xyz = c2h_start(r_hf)
                comment = (f"{xc}/{bas} c2h_saddle E=n/a Eh constructed C2h interchange-saddle start: "
                           f"F..F {C2H_FF:.2f} A, r(H-F) {r_hf:.4f} A from the April acceptor bond, "
                           f"tilt {C2H_TILT_DEG:.0f} deg, exact inversion symmetry")
            d = dist
            print(f"{start:<11}{xc:<6}{bas:<13}{rec['energy'] if start == 'global_min' else float('nan'):>16.9f} "
                  f"{d(xyz[3], xyz[0]):6.3f} {d(xyz[0], xyz[1]):7.3f} {d(xyz[1], xyz[2]):6.3f} "
                  f"{d(xyz[1], xyz[3]):7.3f} {d(xyz[2], xyz[3]):7.3f}  {comment.split(' Eh ', 1)[1][:60]}")
            lines.append(str(NATM))
            lines.append(comment)
            for sym, (x, y, z) in zip(ATOMS, xyz):
                lines.append(f"{sym:<3}{x:16.8f}{y:16.8f}{z:16.8f}")

os.makedirs(os.path.dirname(XYZOUT), exist_ok=True)
open(XYZOUT, "w").write("\n".join(lines) + "\n")
print(f"\nwrote {len(lines) // (NATM + 2)} geometries to {XYZOUT}")
