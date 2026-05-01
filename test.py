import FreeCAD
import Part
import math


# ============================================================
# PARAMETERS (single source of truth)
# ============================================================

class Params:
    """
    Defines a parametric lattice system.

    spacing:
        Period of the lattice (center-to-center distance).

    grid_angle:
        Relative angle between the two lattice families (degrees).

    grid_width:
        Thickness of the SOLID lattice ribs.

    grid_orientation:
        Base orientation of first lattice family (degrees).

    line_width:
        Used only to derive default thickness scaling (optional legacy).
    """

    def __init__(
        self,
        grid_spacing,
        grid_angle,
        grid_width,
        grid_orientation=30,
        line_width=0.4,
    ):
        self.grid_spacing = grid_spacing
        self.grid_angle = grid_angle
        self.grid_width = grid_width

        self.grid_orientation = grid_orientation
        self.line_width=line_width

        self.grid_thickness = self.line_width * 2.1  # optional legacy


# ============================================================
# GRID GENERATOR
# ============================================================

class TiltedGrid:
    def __init__(self, bb, params: Params):
        self.p = params
        self.bb = bb

        self.spacing = params.grid_spacing

        self.xlen = bb.XLength
        self.ylen = bb.YLength
        self.zlen = bb.ZLength

        self.diag_xy = math.sqrt(self.xlen**2 + self.ylen**2)
        self.diag_3d = math.sqrt(self.xlen**2 + self.ylen**2 + self.zlen**2)

        # bounding box for clipping
        self.bbox_solid = Part.makeBox(
            self.xlen,
            self.ylen,
            self.zlen,
            FreeCAD.Vector(bb.XMin, bb.YMin, bb.ZMin),
        )

        # center of bounding box
        self.center = FreeCAD.Vector(
            bb.XMin + self.xlen / 2,
            bb.YMin + self.ylen / 2,
            bb.ZMin + self.zlen / 2,
        )

    # --------------------------------------------------------

    def generate(self, doc, orientation_deg, thickness, tilt_deg=0, phase=0.0):
        """
        Generates one lattice family.

        orientation_deg:
            In-plane rotation of lattice direction.

        thickness:
            Rib thickness (solid or void phase).

        tilt_deg:
            Out-of-plane tilt of lattice.

        phase:
            0.0 = aligned grid
            0.5 = half-period shifted (complementary/negative grid)
        """

        num_slabs = int(self.diag_xy / self.spacing) + 3

        for i in range(-num_slabs, num_slabs):
            offset = (i + phase) * self.spacing

            # --- slab in local coordinates ---
            slab = Part.makeBox(
                thickness,
                self.diag_3d * 2,
                self.diag_3d * 2,
                FreeCAD.Vector(
                    -thickness / 2,
                    -self.diag_3d,
                    -self.diag_3d,
                ),
            )

            # translate along lattice axis
            slab.translate(FreeCAD.Vector(offset, 0, 0))

            # in-plane rotation
            slab.rotate(
                FreeCAD.Vector(0, 0, 0),
                FreeCAD.Vector(0, 0, 1),
                orientation_deg,
            )

            # out-of-plane tilt
            if tilt_deg != 0:
                theta = math.radians(orientation_deg)

                axis = FreeCAD.Vector(
                    -math.sin(theta),
                    math.cos(theta),
                    0,
                )

                slab.rotate(
                    FreeCAD.Vector(0, 0, 0),
                    axis,
                    tilt_deg,
                )

            # move into position
            slab.translate(self.center)

            # clip to bounding box
            slab = slab.common(self.bbox_solid)

            if not slab.isNull():
                obj = doc.addObject(
                    "Part::Feature",
                    f"Lattice_o{orientation_deg:.1f}_t{tilt_deg:.1f}_{i}",
                )
                obj.Shape = slab

        doc.recompute()


# ============================================================
# BUILD FULL STRUCTURE
# ============================================================

def make_grid(doc, bb, p: Params):
    grid = TiltedGrid(bb, p)

    # derived thickness split
    solid = p.grid_thickness
    void = p.grid_spacing -  p.grid_thickness - p.grid_width
    print(f'-----------')
    print(f'grid_spacing {p.grid_spacing}')
    print(f'void {void}')
    print(f'-----------')

    base = p.grid_orientation
    offset = p.grid_angle

    # --------------------------------------------------------
    # FAMILY 1
    # --------------------------------------------------------
    grid.generate(doc, orientation_deg=base, thickness=solid, tilt_deg=30, phase=0.0)
    grid.generate(doc, orientation_deg=base, thickness=void,  tilt_deg=30, phase=0.5)

    # --------------------------------------------------------
    # FAMILY 2 (rotated)
    # --------------------------------------------------------

    grid.generate(doc, orientation_deg=base + offset, thickness=solid, tilt_deg=-30, phase=0.0)
    grid.generate(doc, orientation_deg=base + offset, thickness=void,  tilt_deg=-30, phase=0.5)


# ============================================================
# MAIN
# ============================================================

def see_objects(proj):
    for o in proj.Objects:
        print(o.Name, "|", o.Label, "|", o.TypeId)


if __name__ == "__main__":

    p = Params(
        grid_spacing=30,
        line_width = 0.4,
        grid_angle=0,
        grid_orientation=0,
        grid_width = 5
    )

    project_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    proj = FreeCAD.open(project_path)

    see_objects(proj)

    full_wing = proj.getObject("Pad")
    bb = full_wing.Shape.BoundBox

    make_grid(proj, bb, p)