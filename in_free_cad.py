import FreeCAD
import Part
import math as m


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
        void_scaling=0.9
    ):
        self.grid_spacing = grid_spacing
        self.grid_angle = grid_angle
        self.grid_width = grid_width

        self.grid_orientation = grid_orientation
        self.line_width=line_width
        self.void_scaling = void_scaling

        self.grid_thickness = self.line_width * 2.1  # optional legacy


# ============================================================
# GRID GENERATOR
# ============================================================

class TiltedGrid:
    def __init__(self, doc, bb, shape, params: Params):
        self.p = params
        self.bb = bb
        self.doc = doc

        self.spacing = params.grid_spacing

        self.xlen = bb.XLength
        self.ylen = bb.YLength
        self.zlen = bb.ZLength

        self.diag_xy = m.sqrt(self.xlen**2 + self.ylen**2)
        self.diag_3d = m.sqrt(self.xlen**2 + self.ylen**2 + self.zlen**2)

        # bounding box for clipping
        self.bbox_solid = Part.makeBox(
            self.xlen,
            self.ylen,
            self.zlen,
            FreeCAD.Vector(bb.XMin, bb.YMin, bb.ZMin),
        )

        self.shape = shape.Shape if hasattr(shape, "Shape") else shape

        # center of bounding box
        self.center = FreeCAD.Vector(
            bb.XMin + self.xlen / 2,
            bb.YMin + self.ylen / 2,
            bb.ZMin + self.zlen / 2,
        )

    # --------------------------------------------------------

    def make_grid(self):
        # derived thickness split
        solid = self.p.grid_thickness
        void = self.p.grid_spacing -  self.p.grid_thickness - self.p.grid_width
        print(f'-----------')
        print(f'grid_spacing {self.p.grid_spacing}')
        print(f'void {void}')
        print(f'-----------')

        base = self.p.grid_orientation
        offset = self.p.grid_angle

        shapes = []

        s = self.generate_slabs(orientation_deg=base, thickness=solid, tilt_deg=30, phase=0.0)
        if s: shapes.append(s)
        u1 = (m.cos(base), m.sin(base), 0)
        s = self.generate_slabs(orientation_deg=base, thickness=void, tilt_deg=30, phase=0.5, scale=self.p.void_scaling, u=u1)
        if s: shapes.append(s)

        s = self.generate_slabs(orientation_deg=base + offset, thickness=solid, tilt_deg=-30, phase=0.0)
        if s: shapes.append(s)
        u2 = (m.cos(base + self.p.grid_angle), m.sin(base + self.p.grid_angle), 0)
        s = self.generate_slabs(orientation_deg=base + offset, thickness=void, tilt_deg=-30, phase=0.5, scale=self.p.void_scaling, u=u2)
        if s: shapes.append(s)
        final = Part.makeCompound(shapes)

        obj = self.doc.addObject("Part::Feature", "FullGrid")
        obj.Shape = final

        self.doc.recompute()

    def generate_slabs(self, orientation_deg, thickness, tilt_deg=0, phase=0.0, scale=1, u=None):

        num_slabs = int(self.diag_xy / self.spacing) + 3
        slabs = []

        for i in range(-num_slabs, num_slabs):
            offset = (i + phase) * self.spacing

            slab = Part.makeBox(
                thickness,
                self.ylen * 2,
                self.zlen * 2,
                FreeCAD.Vector(-thickness/2, -self.ylen, -self.zlen),
            )

            slab.translate(FreeCAD.Vector(offset, 0, 0))

            slab.rotate(FreeCAD.Vector(0,0,0), FreeCAD.Vector(0,0,1), orientation_deg)

            if tilt_deg != 0:
                theta = m.radians(orientation_deg)
                axis = FreeCAD.Vector(-m.sin(theta), m.cos(theta), 0)
                slab.rotate(FreeCAD.Vector(0,0,0), axis, tilt_deg)

            slab.translate(self.center)
            clipped_slab = slab.common(self.shape)

            # normalize result to a usable shape

            print(f"\n--- slab {i} ---")
            print("before clip solids:", len(slab.Solids))

            if clipped_slab.isNull():
                continue

            if clipped_slab.ShapeType == "Compound":
                parts = clipped_slab.Solids
                if not parts:
                    continue
                clipped_slab = Part.makeCompound(parts)

            print("after clip valid:", not clipped_slab.isNull())
            print("shape type:", clipped_slab.ShapeType)

            if scale != 1.0:
                clipped_slab = self.scale_perpendicular_to_grid(shape=clipped_slab, factor=scale, u=u)
            slabs.append(clipped_slab)

            print("after scale bbox:", clipped_slab.BoundBox)
            print("center:", safe_center(clipped_slab))

        return Part.makeCompound(slabs)
    
    def scale_perpendicular_to_grid(self, shape, factor, u):
        com = safe_center(shape)

        u = FreeCAD.Vector(u)
        u.normalize()

        # build orthonormal frame once per family
        up = FreeCAD.Vector(0, 0, 1)
        if abs(u.dot(up)) > 0.99:
            up = FreeCAD.Vector(0, 1, 0)

        v = up.cross(u)
        v.normalize()

        w = u.cross(v)
        w.normalize()

        # rotation matrix (local → world)
        R = FreeCAD.Matrix()
        R.A11, R.A12, R.A13 = u.x, v.x, w.x
        R.A21, R.A22, R.A23 = u.y, v.y, w.y
        R.A31, R.A32, R.A33 = u.z, v.z, w.z

        Rinv = R.inverse()

        T1 = FreeCAD.Matrix()
        T1.move(-com)

        # IMPORTANT: no scaling along u (slab direction)
        S = FreeCAD.Matrix()
        S.A11 = 1.0
        S.A22 = factor
        S.A33 = factor

        T2 = FreeCAD.Matrix()
        T2.move(com)

        M = T2.multiply(R.multiply(S.multiply(Rinv.multiply(T1))))

        return shape.transformGeometry(M)
    
def safe_center(shape):
    if hasattr(shape, "CenterOfMass"):
        return shape.CenterOfMass

    # fallback for compounds
    try:
        solids = shape.Solids
        if solids:
            return solids[0].CenterOfMass
    except:
        pass

    raise ValueError("No valid geometry for COM")
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
    full_wing = doc.getObject("Pad")
    lowest_face = min(faces, key=lambda f: f.CenterOfMass.z)
    highest_face = max(faces, key=lambda f: f.CenterOfMass.z)

    # --------------------------------------------------------
    # visualization (fixed naming)
    # --------------------------------------------------------

    viz_low = doc.addObject("Part::Feature", "LowestFaceViz")
    viz_low.Shape = lowest_face

    viz_high = doc.addObject("Part::Feature", "HighestFaceViz")
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

        out = doc.addObject("Part::Feature", "HollowWing")
        out.Shape = hollow_wing

    except Exception as e:
        print("Thickness operation failed:", e)

    doc.recompute()  

def see_objects(doc):
    for o in doc.Objects:
        print(o.Name, "|", o.Label, "|", o.TypeId)


if __name__ == "__main__":

    p = Params(
        grid_spacing=30,
        line_width = 0.4,
        grid_angle=0,
        grid_orientation=0,
        grid_width = 0, 
        void_scaling = 0.6
    )

    docect_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(docect_path)

    see_objects(doc)

    full_wing = doc.getObject("Pad")
    bb = full_wing.Shape.BoundBox
    grid = TiltedGrid(bb=bb, params=p, shape=full_wing, doc=doc)
    grid.make_grid()


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
    # inner = doc.addObject("Part::Feature", "InnerShape")
    # inner.Shape = scaled

    doc.recompute()








    
    # make_grid(doc=doc, bb=bb, p=p, shape=full_wing.Shape)
    # tg = ThicknessGenerator(doc)
    # tg.from_object(obj_name="Pad")