"""Render a planform cross-section of the structured example as an SVG.

Run after main.py --preset example:

    & 'C:\\Program Files\\FreeCAD 1.1\\bin\\python.exe' examples/make_preview.py

Slices the finished STL parallel to the planform (perpendicular to the
thickness axis) and writes the outline plus every internal boundary. That
is the view where the structure is legible: the rib grid runs *through*
the thickness, so a planform cut crosses every rib, whereas an airfoil
section mostly does not -- and stays a single closed loop, because the
bridges deliberately keep the part in one piece.

Pure trimesh/shapely -- no FreeCAD, no matplotlib.
"""

import os
import warnings

import numpy as np
import trimesh

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
STL = os.path.join(HERE, "example_wing_structured.stl")
OUT = os.path.join(HERE, "example_preview.svg")


def rings_of(poly):
    yield np.asarray(poly.exterior.coords)
    for interior in poly.interiors:
        yield np.asarray(interior.coords)


def main():
    if not os.path.exists(STL):
        raise SystemExit(f"{STL} not found -- run: python main.py --preset example")

    mesh = trimesh.load(STL, force='mesh')
    lo, hi = mesh.bounds

    # Scan across the thickness and keep the richest slice rather than a
    # fixed one. How much structure a planform cut reveals varies a lot
    # with height -- near the skin it catches mostly rib edges, nearer the
    # midplane it catches the lightening holes too. Picking by ring count
    # keeps the preview honest without hand-tuning a magic constant.
    best = None
    for frac in np.linspace(0.2, 0.8, 25):
        y = lo[1] + frac * (hi[1] - lo[1])
        section = mesh.section(plane_origin=[0, y, 0], plane_normal=[0, 1, 0])
        if section is None:
            continue
        try:
            planar, _ = section.to_planar()
        except Exception:
            continue
        polys = planar.polygons_full
        if not polys:
            continue
        rings = sum(1 + len(p.interiors) for p in polys)
        if best is None or rings > best[0]:
            best = (rings, y, polys)

    if best is None:
        raise SystemExit("no usable section found")
    _, y, polys = best

    paths, n_rings = [], 0
    for poly in polys:
        for ring in rings_of(poly):
            n_rings += 1
            d = " ".join(f"{'M' if i == 0 else 'L'} {x:.3f} {z:.3f}"
                         for i, (x, z) in enumerate(ring))
            paths.append(d + " Z")

    xs = np.concatenate([np.asarray(p.exterior.coords)[:, 0] for p in polys])
    zs = np.concatenate([np.asarray(p.exterior.coords)[:, 1] for p in polys])
    pad = 4.0
    minx, maxx = xs.min() - pad, xs.max() + pad
    minz, maxz = zs.min() - pad, zs.max() + pad
    w, h = maxx - minx, maxz - minz

    body = "\n".join(f'    <path d="{d}"/>' for d in paths)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg"
     viewBox="{minx:.2f} {minz:.2f} {w:.2f} {h:.2f}"
     width="{w * 3:.0f}" height="{h * 3:.0f}">
  <title>ContinuousPath example wing - planform section</title>
  <g fill="#7ee7a8" fill-rule="evenodd" stroke="#123" stroke-width="0.35"
     stroke-linejoin="round">
{body}
  </g>
</svg>
'''
    with open(OUT, 'w', encoding='utf-8') as fh:
        fh.write(svg)

    print(f"section at y={y:.2f}: {len(polys)} polygon(s), {n_rings} ring(s)")
    print(f"wrote {OUT} ({os.path.getsize(OUT) / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
