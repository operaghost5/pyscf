# Why the original CNEO-DFT runs used DFT grid level 3

Note for the CNEO-DFT hydrogen-bonded dimer project, 21 September 2026.
It answers the question why the electronic DFT grids appeared to be reset
from level 9 to level 3 in the original CNEO-DFT calculations (February and
April 2026). File and line references are to `pyscf/neo/ks.py` and
`pyscf/dft/gen_grid.py` in the theorychemyang/pyscf fork as checked out for
this project; the patched copy that ran on the cluster is quoted from the
excerpts shared earlier.

## Short answer

Level 9 was never applied to the electronic component, so nothing was reset.
The electronic Kohn–Sham object inside each NEO calculation kept PySCF's
default grid level 3 from the first SCF cycle to the last Hessian, for three
reasons that reinforced one another:

1. the input files set no grid at all;
2. the `grids.level = 9` lines added to `neo/ks.py` sat in code that only
   runs when an electron–proton correlation (EPC) functional is present, and
   the calculations used none;
3. the `grids` attribute of the NEO object, the one the debug print reported
   as "set to None", is the EPC screening grid, not the electronic grid.

The conventional-DFT runs set `mf.grids.level = 9` directly on a plain
`dft.RKS` object, which is why they ran at level 9 while the CNEO-DFT runs of
the same systems ran at level 3.

## 1. Where PySCF reads the grid in a NEO calculation

`neo.CDFT(mol, xc=...)` builds one Kohn–Sham object per component: the
electrons in `mf.components['e']` (an ordinary `pyscf.dft.rks.RKS`, or `UKS`
for open shells) and one Hartree–Fock-like object per quantum nucleus. The
electronic exchange–correlation integration uses the electronic component's
own `Grids` object, `mf.components['e'].grids`. Its `level` defaults to 3
(`gen_grid.py:574`, config key `dft_gen_grid_Grids_level`), and the grid
points are generated lazily, from whatever `level` or `atom_grid` is set at
that moment, the first time the electronic `get_veff` runs.

The original drivers never touched that object. The calculation section of
every original input reads

```
mol = neo.Mole()
mol.build(atom=gm_start_geom_hydb,
          basis=bas_two, nuc_basis='pb4d', charge = 0, unit='Bohr', quantum_nuc=[0])

mf = neo.CDFT(mol, xc=func_one)
mf.verbose = 10
```

with no grid line, so the electronic component ran at level 3. The NEO
driver does not dump the component's grid flags, so the logs carry no grid
printout; the level-3 conclusion follows from the input and is consistent
with the size of the translational sum-rule violations analysed in the
methods document (for example 0.07 to 0.4 Eh/Bohr² on chlorine, the
signature of 80 radial shells).

## 2. The patched lines sat in EPC-only code

Two places in the cluster copy of `neo/ks.py` were patched. Both are shown
as they were shared, with the early return that made the patch unreachable
marked.

**(a) `InteractionCorrelation.get_vint`**, the per-pair interaction routine.

```
vj = super().get_vint(dm, *args, **kwargs)
# For nuclear initial guess, use Coulomb only
if coulomb_only or \
        not (self.mf1_type in dm and self.mf2_type in dm and self._need_epc()):
    return vj                                   # <- taken on every call when epc is None
...
ao_loc_n = mol_n.ao_loc_nr()

mf_e.grids.level = 9 #GBT added this line
grids_e = mf_e.grids
grids_changed = (self._elec_grids_hash != _hash_grids(grids_e))
```

`_need_epc()` returns `False` when the interaction carries no EPC functional,
so the method returns the Coulomb matrices before it reaches the patched
line (current fork: return at `ks.py:175–177`, grid line at `ks.py:202`).
The logs show `GBT ... calling get_vint() from neo.ks.InteractionCorrelation`
on every SCF cycle, so the method was entered each time; the return at the
EPC check is what stopped the rest of it from running.

**(b) The fast EPC path**, `KS.get_vint_fast` in the cluster copy
(`_get_epc_vmat` in the current fork).

```
vj = super().get_vint(mol, dm, no_epc=True)
if self.epc is None:
    return vj                                   # <- taken on every call when epc is None

mf_e = self.components['e']
grids_e = mf_e.grids
print('GBT ... adding line below which modifies mf_e.grids.level')
mf_e.grids.level = 9 #GBT added this line
grids_e = mf_e.grids #GBT added this line too
print('GBT ... added line above which modifies mf_e.grids.level')
grids_changed = (self._elec_grids_hash != _hash_grids(grids_e))
print('GBT added this line which prints Boolean of grids_changed: ', grids_changed)
```

The three `print` statements around the patched line appear in none of the
logs of any system (zero occurrences across every log set shared in this
project), which confirms directly that this code was never reached. The
current fork keeps the same order: `if mf.epc is None: return no_epc()` at
`ks.py:284` precedes `grids_e = mf_e.grids` at `ks.py:288`.

Even with an EPC functional the placement would have been fragile. Both
routines run inside the SCF, after the electronic grid has already been
generated for the current geometry, and changing `level` on an existing
`Grids` object does not regenerate its points; the new level would have
taken effect only after the next `reset`, that is from the second geometry
step onward.

## 3. `mf.grids` on the NEO object is the EPC screening grid

Most of the CNEO-DFT logs open with the debug line

```
GBT ... neo.KS class initialization complete ... self.grids set to None and self._elec_grids_hash set to None
```

That attribute is not the electronic grid. It is set to `None` in
`InteractionCorrelation.__init__` (`ks.py:141–142`), in `KS.__init__`
(`ks.py:483–484`) and again at the end of `KS.reset` (`ks.py:658–659`), which
the geometry optimizer calls with the new geometry at every step. When an
EPC functional is present, `get_vint` fills it with a screened copy of the
electronic grid restricted to the blocks where protonic basis functions are
non-zero (`ks.py:209–217`), and uses that copy only for the EPC integrals.
It never enters the electronic exchange–correlation integration. Two
consequences:

- "self.grids set to None" in the log is not evidence that the electronic
  grid was reset; the electronic `Grids` object lives on the component and
  was never touched.
- a level assigned to `mf.grids` on the NEO object (after creating a `Grids`
  object there) would be discarded at the first optimizer step and would in
  any case not be read by the electronic integration.

## 4. The setting that works

The grid must be set on the electronic component, before the optimization
starts:

```
mf = neo.CDFT(mol, xc=XC)
mf.components['e'].grids.level = 5
mf.components['e'].grids.atom_grid = {'Br': (400, 770)}   # dense radial grid on the halogen
```

This survives the per-step `reset`. `KS.reset` keeps the existing component
objects and calls their own `reset(comp)` (`ks.py:626`), and `Grids.reset`
(`gen_grid.py:687–697`) clears only the generated data (`coords`, `weights`,
`non0tab`, `screen_index`, `atm_idx`, `quadrature_weights`) while keeping
`level`, `atom_grid`, `prune` and `radi_method`. The grid is therefore
regenerated at every geometry with the requested settings. All the new
CNEO-DFT drivers in `nqe_dimer_inputs` set the grid this way and print a
grid report after the first SCF (level, radial scheme, partition, pruning,
radial × angular points per element, total points), so the applied grid is
visible in every log; the conventional-DFT drivers set the same values on
`mf.grids` of the `dft.RKS` object and print the same report.

## 5. Grid sizes for reference

Radial shells × Lebedev angular points per atom before NWChem pruning, from
the `RAD_GRIDS` and `ANG_ORDER` tables (`gen_grid.py:753–780`).

| Element | level 3 (original CNEO-DFT) | level 5 (new inputs) | level 9 (original conventional DFT) | dense halogen grid (new inputs) |
|---|---|---|---|---|
| H | 50 × 302 | 70 × 590 | 200 × 1454 | – |
| C, N, F | 75 × 302 | 105 × 770 | 200 × 1454 | F: 400 × 770 in the HF sets and the `f400` variants |
| Cl | 80 × 434 | 110 × 770 | 200 × 1454 | 400 × 770 |
| Br | 90 × 434 | 120 × 770 | 200 × 1454 | 400 × 770 |

The level-9 entry is the setting seen in the conventional-DFT inputs shared
for this project; the HCN···HBr conventional-DFT input was not among them
and is assumed to match.
