"""Generate the example wing shipped with this repository.

Produces a small, tapered, swept wing as a proper BREP solid by lofting
NACA 4-digit airfoil sections. Run with FreeCAD's bundled Python:

    & 'C:\\Program Files\\FreeCAD 1.1\\bin\\python.exe' examples/make_example_wing.py

The geometry is generated from the published NACA 4-digit formula, so the
example is original to this project and carries no third-party model
licence. It is also deliberately a *lofted* solid rather than a
mesh-derived one: it has a handful of real surfaces instead of one planar
face per triangle, which is what a well-formed input to this tool looks
like (see the "Input requirements" section of the README).

The wing stands along Z so it slices the way the presets expect:
    X = chord, Y = thickness, Z = span
"""

import os
import sys

import FreeCAD
import Part

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "example_wing.FCStd")
OBJ_NAME = "ExampleWing"


def naca4(thickness_pct=12.0, camber_pct=2.0, camber_pos=0.4, n=32):
    """NACA 4-digit section as a closed list of (x, y) in unit chord.

    Standard published formulation: a thickness distribution yt applied
    either side of a two-segment parabolic camber line. Points are
    cosine-spaced so the leading edge, where curvature is highest, gets
    the most resolution for a given point count.
    """
    import math

    t = thickness_pct / 100.0
    m = camber_pct / 100.0
    p = camber_pos

    # cosine spacing: dense at the leading edge, sparse at the trailing
    xs = [(1 - math.cos(math.pi * i / (n - 1))) / 2.0 for i in range(n)]

    upper, lower = [], []
    for x in xs:
        yt = 5 * t * (0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                      + 0.2843 * x ** 3 - 0.1015 * x ** 4)
        if m > 0 and 0 < x < 1:
            if x < p:
                yc = m / p ** 2 * (2 * p * x - x ** 2)
                dyc = 2 * m / p ** 2 * (p - x)
            else:
                yc = m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2)
                dyc = 2 * m / (1 - p) ** 2 * (p - x)
            th = math.atan(dyc)
            upper.append((x - yt * math.sin(th), yc + yt * math.cos(th)))
            lower.append((x + yt * math.sin(th), yc - yt * math.cos(th)))
        else:
            upper.append((x, yt))
            lower.append((x, -yt))

    # walk the upper surface trailing->leading, then the lower back again
    pts = list(reversed(upper)) + lower[1:]

    # The NACA thickness function leaves a small gap at x=1; close it so
    # the lofted solid is watertight rather than needing a repair pass.
    if pts[0] != pts[-1]:
        pts.append(pts[0])
    return pts


def section_wire(chord, span_pos, sweep, twist_deg=0.0):
    """One airfoil section placed at its spanwise station."""
    import math

    pts2d = naca4()
    tw = math.radians(twist_deg)
    verts = []
    for x, y in pts2d:
        xc, yc = x * chord, y * chord
        # twist about the quarter-chord
        qx = xc - 0.25 * chord
        rx = qx * math.cos(tw) - yc * math.sin(tw) + 0.25 * chord
        ry = qx * math.sin(tw) + yc * math.cos(tw)
        verts.append(FreeCAD.Vector(rx + sweep, ry, span_pos))
    return Part.makePolygon(verts)


def build_wing(root_chord=100.0, tip_chord=55.0, span=200.0,
               sweep=35.0, twist=-3.0):
    """Loft three sections into a closed solid.

    A mid station is included so the taper is not perfectly linear --
    that gives the rib grid something more interesting to cut than a
    plain wedge.
    """
    sections = [
        section_wire(root_chord, 0.0, 0.0, 0.0),
        section_wire(root_chord * 0.72 + tip_chord * 0.28,
                     span * 0.55, sweep * 0.5, twist * 0.5),
        section_wire(tip_chord, span, sweep, twist),
    ]
    # solid=True, ruled=False -- a smooth lofted surface rather than one
    # planar strip per point pair. Both are valid, but the smooth loft
    # produces 65 faces instead of 192 and a 286KB STEP instead of 672KB,
    # and it is the better example: a handful of real surfaces is exactly
    # the well-formed input this tool wants, and what makes the BREP
    # fallback path cheap rather than pathological.
    solid = Part.makeLoft(sections, True, False)
    if not solid.isValid():
        solid = solid.removeSplitter()
    return solid


def main():
    wing = build_wing()
    bb = wing.BoundBox
    print(f"faces={len(wing.Faces)} solids={len(wing.Solids)} "
          f"shells={len(wing.Shells)}")
    print(f"valid={wing.isValid()} closed={wing.isClosed()} "
          f"volume={wing.Volume:.1f} mm^3")
    print(f"bbox: X={bb.XLength:.1f} Y={bb.YLength:.1f} Z={bb.ZLength:.1f}")

    if not (wing.isClosed() and wing.Volume > 0 and len(wing.Shells) == 1):
        print("ERROR: generated wing is not a closed single-shell solid")
        return 1

    # Saved as a FreeCAD document rather than STEP on purpose. STEP has no
    # single authoritative unit -- FreeCAD's importer applies a preference
    # that varies per machine, and a round-trip here came back scaled by
    # 25.4 (mm read as inches), turning a 200mm wing into a 5m one. An
    # example that silently depends on a local setting is worse than no
    # example, and .FCStd carries its units unambiguously.
    doc = FreeCAD.newDocument("ExampleWing")
    feature = doc.addObject("Part::Feature", OBJ_NAME)
    feature.Shape = wing
    doc.recompute()
    doc.saveAs(OUT)

    print(f"wrote {OUT} ({os.path.getsize(OUT) / 1024:.0f} KB), "
          f"object '{OBJ_NAME}'")

    # Reload and confirm the round-trip preserves scale -- the check that
    # would have caught the STEP unit problem immediately.
    reread = FreeCAD.open(OUT).getObject(OBJ_NAME).Shape
    if abs(reread.Volume - wing.Volume) / wing.Volume > 1e-6:
        print(f"ERROR: round-trip changed the volume "
              f"({wing.Volume:.1f} -> {reread.Volume:.1f})")
        return 1
    rb = reread.BoundBox
    print(f"round-trip OK: volume={reread.Volume:.1f} "
          f"bbox={rb.XLength:.1f} x {rb.YLength:.1f} x {rb.ZLength:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
