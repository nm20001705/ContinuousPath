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
    def __init__(self, bb, shape, params: Params):
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

        self.shape = shape

        # center of bounding box
        self.center = FreeCAD.Vector(
            bb.XMin + self.xlen / 2,
            bb.YMin + self.ylen / 2,
            bb.ZMin + self.zlen / 2,
        )

    # --------------------------------------------------------

    def generate(self, doc, orientation_deg, thickness, tilt_deg=0, phase=0.0):

        num_slabs = int(self.diag_xy / self.spacing) + 3
        clipped = []

        # Precompute rotation
        rot_main = FreeCAD.Rotation(FreeCAD.Vector(0,0,1), orientation_deg)

        # Tilt axis (same as your code)
        if tilt_deg != 0:
            theta = math.radians(orientation_deg)
            tilt_axis = FreeCAD.Vector(-math.sin(theta), math.cos(theta), 0)
            rot_tilt = FreeCAD.Rotation(tilt_axis, tilt_deg)
        else:
            rot_tilt = None

        for i in range(-num_slabs, num_slabs):
            offset = (i + phase) * self.spacing

            # ------------------------------------------------------------
            # 1. Create a reasonably sized slab (NOT gigantic)
            # ------------------------------------------------------------
            slab = Part.makeBox(
                thickness,
                self.ylen * 2,
                self.zlen * 2,
                FreeCAD.Vector(-thickness/2, -self.ylen, -self.zlen)
            )

            # ------------------------------------------------------------
            # 2. Apply main rotation
            # ------------------------------------------------------------
            slab = slab.rotate(FreeCAD.Vector(0,0,0), FreeCAD.Vector(0,0,1), orientation_deg)

            # ------------------------------------------------------------
            # 3. Apply tilt (if any)
            # ------------------------------------------------------------
            if rot_tilt:
                slab = slab.rotate(FreeCAD.Vector(0,0,0), tilt_axis, tilt_deg)

            # ------------------------------------------------------------
            # 4. Translate slab along X (your conceptual logic)
            # ------------------------------------------------------------
            slab.translate(FreeCAD.Vector(offset, 0, 0))

            # ------------------------------------------------------------
            # 5. Move slab to bounding box center
            # ------------------------------------------------------------
            slab.translate(self.center)

            # ------------------------------------------------------------
            # 6. Clip slab individually (CRITICAL FIX)
            # ------------------------------------------------------------
            clipped_slab = slab.common(self.bbox_solid)

            if clipped_slab.isNull():
                continue
            if not clipped_slab.isValid():
                continue

            clipped.append(clipped_slab)

        # ------------------------------------------------------------
        # 7. Fuse only valid clipped slabs
        # ------------------------------------------------------------
        if not clipped:
            return None

        final = Part.makeCompound(clipped)

        if not final.isValid():
            print("WARNING: final compound invalid")
            return None

        return final




# ============================================================
# BUILD FULL STRUCTURE
# ============================================================

def make_grid(doc, bb, shape, p: Params):
    grid = TiltedGrid(bb=bb, params=p, shape=shape)

    # derived thickness split
    solid = p.grid_thickness
    void = p.grid_spacing -  p.grid_thickness - p.grid_width
    print(f'-----------')
    print(f'grid_spacing {p.grid_spacing}')
    print(f'void {void}')
    print(f'-----------')

    base = p.grid_orientation
    offset = p.grid_angle

    shapes = []

    s = grid.generate(doc, orientation_deg=base, thickness=solid, tilt_deg=30, phase=0.0)
    if s: shapes.append(s)

    s = grid.generate(doc, orientation_deg=base, thickness=void, tilt_deg=30, phase=0.5)
    if s: shapes.append(s)

    s = grid.generate(doc, orientation_deg=base + offset, thickness=solid, tilt_deg=-30, phase=0.0)
    if s: shapes.append(s)

    s = grid.generate(doc, orientation_deg=base + offset, thickness=void, tilt_deg=-30, phase=0.5)
    if s: shapes.append(s)
    final = Part.makeCompound(shapes)

    obj = doc.addObject("Part::Feature", "FullGrid")
    obj.Shape = final

    doc.recompute()

# ============================================================
# FACE ANALYSIS
# ============================================================

class FaceAnalyzer:
    """
    Utility class for geometric face inspection.
    """

    @staticmethod
    def centroid(face):
        return face.CenterOfMass

    @staticmethod
    def normal(face):
        n = face.normalAt(0.5, 0.5)
        n.normalize()
        return n


# ============================================================
# THICKNESS BUILDER
# ============================================================

class ThicknessBuilder:
    """
    Builds a solid from one or two selected faces.

    Rules:
    - 1 face:
        → extrude along face normal
    - 2 faces:
        → bottom = lowest centroid Z
        → top = highest centroid Z
        → build enclosed solid
    """

    def __init__(self, shape):
        self.shape = shape
        self.faces = shape.Faces

        if len(self.faces) == 0:
            raise ValueError("Shape has no faces")

    # --------------------------------------------------------

    def build(self, thickness=None):
        """
        Returns a solid based on face configuration.
        """

        if len(self.faces) == 1:
            return self._from_single_face(thickness)

        if len(self.faces) >= 2:
            return self._from_two_faces()

        raise ValueError("Unsupported geometry")

    # --------------------------------------------------------

    def _from_single_face(self, thickness):
        face = self.faces[0]

        n = FaceAnalyzer.normal(face)

        if thickness is None:
            thickness = 1.0

        solid = face.extrude(n.multiply(thickness))
        return solid

    # --------------------------------------------------------

    def _from_two_faces(self):
        # sort faces by centroid Z
        sorted_faces = sorted(
            self.faces,
            key=lambda f: FaceAnalyzer.centroid(f).z
        )

        bottom = sorted_faces[0]
        top = sorted_faces[-1]

        # attempt shell construction
        shell = Part.makeShell([bottom, top])

        # convert to solid
        solid = Part.Solid(shell)

        return solid


# ============================================================
# HIGH-LEVEL API WRAPPER
# ============================================================

class ThicknessGenerator:
    """
    User-facing interface.
    """

    def __init__(self, doc):
        self.doc = doc

    def from_object(self, obj_name, out_name="ThickSolid", thickness=None):
        obj = self.doc.getObject(obj_name)

        builder = ThicknessBuilder(obj.Shape)
        solid = builder.build(thickness=thickness)

        out = self.doc.addObject("Part::Feature", out_name)
        out.Shape = solid

        self.doc.recompute()
        return out

# ============================================================
# MAIN
# ============================================================

def make_thickness():
    full_wing = proj.getObject("Pad")
    lowest_face = min(faces, key=lambda f: f.CenterOfMass.z)
    highest_face = max(faces, key=lambda f: f.CenterOfMass.z)

    # --------------------------------------------------------
    # visualization (fixed naming)
    # --------------------------------------------------------

    viz_low = proj.addObject("Part::Feature", "LowestFaceViz")
    viz_low.Shape = lowest_face

    viz_high = proj.addObject("Part::Feature", "HighestFaceViz")
    viz_high.Shape = highest_face

    # --------------------------------------------------------
    # thickness operation
    # --------------------------------------------------------

    try:
        hollow_wing = full_wing.Shape.makeThickness(
            [lowest_face, highest_face],
            0.2,
            0.001
        )

        out = proj.addObject("Part::Feature", "HollowWing")
        out.Shape = hollow_wing

    except Exception as e:
        print("Thickness operation failed:", e)

    proj.recompute()  

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

    # center = full_wing.Shape.CenterOfMass
    # factor = 0.9  # 90% size
    # mat = FreeCAD.Matrix()

    # mat.A11 = factor
    # mat.A22 = factor
    # mat.A33 = factor

    # scaled = full_wing.Shape.copy()
    # scaled = scaled.transformGeometry(mat)

    # # --------------------------------------------------------
    # # move back to center correctly
    # # --------------------------------------------------------
    # shift = FreeCAD.Vector(
    #     center.x * (1 - factor),
    #     center.y * (1 - factor),
    #     center.z * (1 - factor),
    # )

    # scaled.translate(-shift)

    # # --------------------------------------------------------
    # # show result
    # # --------------------------------------------------------
    # inner = proj.addObject("Part::Feature", "InnerShape")
    # inner.Shape = scaled

    proj.recompute()








    
    make_grid(doc=proj, bb=bb, p=p, shape=full_wing.Shape)
    tg = ThicknessGenerator(proj)

    tg.from_object(obj_name="Pad")