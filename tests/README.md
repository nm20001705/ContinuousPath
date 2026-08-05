# Tests

    python -m pip install -r requirements-dev.txt
    python -m pytest

**These run under any Python — FreeCAD is not required.** That is
deliberate, and it is also a real constraint on what they cover.

## What is covered

The pure-geometry core, which imports only numpy / trimesh / shapely:

| module | what is pinned |
|---|---|
| `rib_slice_core` | bridge & hole width policies, the `(u,v)` basis, edge clustering, interval collection |
| `mesh_simplify_utils` | zero-volume sliver removal |

Each test targets a defect that actually occurred, and the docstrings say
which. The recurring pattern was a bug in a few lines of arithmetic that
could only be *observed* after a multi-minute FreeCAD run, so the symptom
appeared far from the cause:

- **Bridges skipped instead of clamped** where a section was thinner than
  `bridge_height`, breaking the bridge column exactly where the wall is
  weakest.
- **Hole margin dropping the `thickness/2` term**, so holes cleared the
  neighbouring rib's centre plane rather than its extruded volume.
- **The hole profile mapped onto the raw span** instead of the
  margin-inset one, shifting the taper.
- **`t` stored absolutely rather than relative to `t_rib`**, which made
  bridges and holes drift sideways as they climbed.
- **Sliver cleanup opening a watertight mesh** — a leaky cut tool sews
  into a shell, `Part.makeSolid()` fails, and the wing cut then returns a
  Null shape or runs for hours.

The tests were verified by mutation: each was re-run against a
reintroduction of the original bug and confirmed to fail. Two cases
needed sharper fixtures before they discriminated at all —
- sliver debris had to be given *more faces than a small solid piece*,
  otherwise filtering on face count and on volume agree by accident;
- the drift tests needed a **tilted** rib, because with an upright rib
  `t_rib` is zero and storing raw `t` is indistinguishable from storing
  the offset.

Both are noted in the tests, since the naive fixture silently proves
nothing.

## What is NOT covered

Everything that needs FreeCAD, real geometry, or the boolean backend:

- input loading, tessellation, and the mesh→BREP conversion
  (`slab_utils`, `assembly_utils_freecad`)
- rib segmentation against a real wing (`build_rib_segments_analytical`)
- the final assembly booleans and STL export
- `viz_utils` (writes into a FreeCAD document)

These are exercised by running the pipeline end to end:

    & 'C:\Program Files\FreeCAD 1.1\bin\python.exe' main.py --preset wingR1

The numbers to compare against a known-good run are the rib-segment
count, the bridge/hole solid face counts, and the final mesh face count
and volume. A pure refactor must leave all of them byte-identical.
