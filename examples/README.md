# Example

A small tapered, swept wing you can run the pipeline on straight from a
fresh clone.

```powershell
& 'C:\Program Files\FreeCAD 1.1\bin\python.exe' main.py --preset example
```

Roughly 10 seconds. It writes `example_wing_structured.stl` next to the
input.

## What you should see

```
Input 'ExampleWing': 65 faces, valid=True, closed=True, solids=1, shells=1,
                     volume=123588.8
Created 77 rib segments (analytical grid).
Bridge solid: 11328 faces, watertight=True
Hole solid:    9308 faces, watertight=True
  rib - bridges: 12348 faces, watertight=True
  rib - holes:   20852 faces, watertight=True
Final (trimesh): 38378 faces, watertight=True, volume=120722.2
  material removed by ribs: 1653.0 mm^3 (1.35% of the wing)
```

`watertight=True` at every stage is the thing to check. If the input
report says anything other than `closed=True` with one shell, stop —
see the README's *Input requirements*.

## Looking at the result

```powershell
python examples/make_preview.py
```

Writes `example_preview.svg`, a planform section through the finished
part.

**Slice it in the right plane or you will think nothing happened.** The
rib grid runs *through* the thickness, so:

| section plane | what it shows |
|---|---|
| planform (perpendicular to thickness) | the full grid — 6 outlines, 29 boundaries |
| airfoil (perpendicular to span) | a single closed loop |

The airfoil section staying one loop is not a failure — it is the
bridges working. They deliberately keep the part in one piece so the
nozzle never has to retract. Only 1.35% of the material is gone; the
structure is the *slit pattern*, not a hollowed interior.

## Regenerating the input

`example_wing.FCStd` is committed, but it is generated, not hand-modelled:

```powershell
& 'C:\Program Files\FreeCAD 1.1\bin\python.exe' examples/make_example_wing.py
```

It lofts NACA 4-digit sections (root chord 100 mm, tip 55 mm, 200 mm
span, 35 mm sweep, −3° washout) into a closed 65-face solid. The airfoil
comes from the published NACA formula, so the example is original to this
project and carries no third-party model licence.

Two deliberate choices, both learned the hard way:

- **Saved as `.FCStd`, not STEP.** STEP has no single authoritative unit.
  FreeCAD's importer applies a per-machine preference, and a round-trip
  here came back scaled by 25.4 — turning the 200 mm wing into a 5 m one
  and the run into a multi-thousand-slice crawl. The generator now
  re-opens what it wrote and fails if the volume moved.
- **A smooth loft, not a ruled one.** 65 real surfaces instead of 192
  planar strips, and 22 KB instead of 672 KB. It is also the better
  example: a handful of genuine surfaces is what well-formed input to
  this tool looks like, as opposed to the mesh-derived solids that make
  the BREP path pathological.

## A note on previews

The `example` preset turns the `vis_*` flags off. The pipeline writes its
preview meshes into the input document and saves it, which took
`example_wing.FCStd` from 22 KB to 741 KB and left it dirty in git.

To inspect the intermediate geometry in the FreeCAD GUI, turn the flags
on in `config.toml`, run, then reset the document:

```powershell
git checkout examples/example_wing.FCStd
```

Even with previews off, FreeCAD rewrites a timestamp on save, so the file
may show as modified by a few bytes after a run. That is harmless.
