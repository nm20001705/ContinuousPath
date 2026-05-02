import FreeCAD
import Part
import math as m


# ============================================================
# PARAMETERS
# ============================================================

class Params:
    """
    Defines a parametric lattice system.

    grid_spacing:
        Period of the lattice (center-to-center distance).

    grid_angle:
        Relative angle between the two lattice families (degrees).

    grid_width:
        Thickness of the SOLID lattice ribs.

    grid_orientation:
        Base orientation of first lattice family (degrees).

    rib_clearance:
        Absolute gap (in model units) left between a void slab and its
        neighboring solid ribs on each side. The void slab is made
        narrower by 2 * rib_clearance before clipping, so the gap is
        constant regardless of how the clip shape cuts the slab.
        Must satisfy: 2 * rib_clearance < void_thickness.
    """

    def __init__(
        self,
        grid_spacing: float,
        grid_angle: float,
        grid_width: float,
        grid_orientation: float = 30.0,
        rib_clearance: float = 1.0,
    ):
        if grid_width >= grid_spacing:
            raise ValueError("grid_width must be less than grid_spacing")
        if rib_clearance < 0:
            raise ValueError("rib_clearance must be >= 0")

        self.grid_spacing = grid_spacing
        self.grid_angle = grid_angle
        self.grid_width = grid_width
        self.grid_orientation = grid_orientation
        self.rib_clearance = rib_clearance

    @property
    def void_thickness(self) -> float:
        """Full width of the void slot between two solid ribs."""
        return self.grid_spacing - self.grid_width

    @property
    def void_slab_thickness(self) -> float:
        """
        Actual thickness of the void slab after subtracting clearance on
        both sides. This is what gets passed to _generate_slabs.
        """
        t = self.void_thickness - 2.0 * self.rib_clearance
        if t <= 0:
            raise ValueError(
                f"rib_clearance ({self.rib_clearance}) is too large: "
                f"void_slab_thickness would be {t:.3f}. "
                f"Reduce rib_clearance or increase grid_spacing."
            )
        return t


# ============================================================
# HELPERS
# ============================================================

def safe_center(shape: Part.Shape) -> FreeCAD.Vector:
    """Return center of mass, with fallback for compounds."""
    if hasattr(shape, "CenterOfMass"):
        return shape.CenterOfMass
    try:
        solids = shape.Solids
        if solids:
            return solids[0].CenterOfMass
    except Exception:
        pass
    raise ValueError("Cannot compute center of mass for shape")


def is_valid_compound(shape) -> bool:
    """Return True if shape is non-null and contains at least one solid."""
    if shape is None or shape.isNull():
        return False
    try:
        return len(shape.Solids) > 0
    except Exception:
        return False


def fuse_compounds(a: Part.Shape, b: Part.Shape) -> Part.Shape:
    """
    Fuse two compounds safely.
    Falls back to a plain compound if the boolean fuse fails,
    which is acceptable for use as a cutting tool.
    """
    try:
        return a.fuse(b)
    except Exception as e:
        print(f"[warn] fuse failed, falling back to compound: {e}")
        return Part.makeCompound(a.Solids + b.Solids)


# ============================================================
# GRID GENERATOR
# ============================================================

class TiltedGrid:
    """
    Generates a two-family tilted slab lattice clipped to an arbitrary shape.

    Each family produces:
      - solid ribs  (grid_width thick)
      - void slabs  (void_slab_thickness wide = void_thickness - 2*rib_clearance,
                     sized before clipping so the gap to neighboring ribs is
                     always approximately rib_clearance in the stacking direction)

    The void compounds from both families are then used to cut the solid ribs.
    """

    def __init__(self, doc, bb, shape, params: Params):
        self.p = params
        self.bb = bb
        self.doc = doc

        self.xlen = bb.XLength
        self.ylen = bb.YLength
        self.zlen = bb.ZLength

        self.diag_xy = m.sqrt(self.xlen ** 2 + self.ylen ** 2)

        self.shape = shape.Shape if hasattr(shape, "Shape") else shape

        self.center = FreeCAD.Vector(
            bb.XMin + self.xlen / 2,
            bb.YMin + self.ylen / 2,
            bb.ZMin + self.zlen / 2,
        )

        # Shrunk clip shape for voids — computed once, reused for all families.
        self.inset_shape = self._make_inset_shape()

    # --------------------------------------------------------
    # Public entry point
    # --------------------------------------------------------

    def make_grid(self):
        base = self.p.grid_orientation
        offset = self.p.grid_angle

        # --- Family 1 ---
        ribs1 = self._make_ribs(orientation_deg=base, tilt_deg=30)
        voids1 = self._make_voids(orientation_deg=base, tilt_deg=30)

        # --- Family 2 ---
        ribs2 = self._make_ribs(orientation_deg=base + offset, tilt_deg=-30)
        voids2 = self._make_voids(orientation_deg=base + offset, tilt_deg=-30)

        # --- Cut ribs by all voids ---
        all_voids = self._combine_voids(voids1, voids2)
        cut_ribs1 = self._cut_safe(ribs1, all_voids)
        cut_ribs2 = self._cut_safe(ribs2, all_voids)

        parts = [s for s in [cut_ribs1, cut_ribs2] if is_valid_compound(s)]
        if not parts:
            print("[warn] make_grid: no valid geometry produced")
            return

        final = Part.makeCompound(parts)
        obj = self.doc.addObject("Part::Feature", "FullGrid")
        obj.Shape = final
        self.doc.recompute()

    # --------------------------------------------------------
    # Private helpers
    # --------------------------------------------------------

    def _make_inset_shape(self) -> Part.Shape:
        """
        Return a uniformly scaled-down copy of self.shape, shrunk by
        approximately rib_clearance on all sides.

        Scale factor: 1 - rib_clearance / characteristic_radius
        where characteristic_radius = half the 3D bounding box diagonal.

        This is a uniform scale from the shape's center of mass, so the
        inset distance is only approximate — it will be exact at points
        that are exactly characteristic_radius from the center, and will
        vary elsewhere. For a convex, roughly symmetric shape like a wing
        this is a good enough approximation.
        """
        com = self.shape.CenterOfMass
        char_radius = m.sqrt(self.xlen**2 + self.ylen**2 + self.zlen**2) / 2.0
        scale = 1.0 - self.p.rib_clearance / char_radius

        if scale <= 0:
            raise ValueError(
                f"rib_clearance ({self.p.rib_clearance}) is too large relative "
                f"to the shape (characteristic_radius={char_radius:.1f}). "
                f"The inset shape would vanish."
            )

        # Translate to origin, scale, translate back
        T1 = FreeCAD.Matrix()
        T1.move(-com)

        S = FreeCAD.Matrix()
        S.A11 = scale
        S.A22 = scale
        S.A33 = scale

        T2 = FreeCAD.Matrix()
        T2.move(com)

        M = T2.multiply(S.multiply(T1))
        return self.shape.transformGeometry(M)

    def _make_ribs(self, orientation_deg: float, tilt_deg: float) -> Part.Shape:
        return self._generate_slabs(
            orientation_deg=orientation_deg,
            thickness=self.p.grid_width,
            tilt_deg=tilt_deg,
            phase=0.0,
            clip_shape=self.shape,
        )

    def _make_voids(self, orientation_deg: float, tilt_deg: float) -> Part.Shape:
        return self._generate_slabs(
            orientation_deg=orientation_deg,
            thickness=self.p.void_slab_thickness,
            tilt_deg=tilt_deg,
            phase=0.5,
            clip_shape=self.inset_shape,
        )

    def _combine_voids(self, voids1, voids2) -> Part.Shape | None:
        valid1 = is_valid_compound(voids1)
        valid2 = is_valid_compound(voids2)
        if valid1 and valid2:
            return fuse_compounds(voids1, voids2)
        if valid1:
            return voids1
        if valid2:
            return voids2
        return None

    def _cut_safe(self, ribs: Part.Shape, voids: Part.Shape) -> Part.Shape | None:
        if not is_valid_compound(ribs):
            return None
        if not is_valid_compound(voids):
            return ribs
        try:
            return ribs.cut(voids)
        except Exception as e:
            print(f"[warn] cut failed: {e}")
            return ribs

    def _generate_slabs(
        self,
        orientation_deg: float,
        thickness: float,
        tilt_deg: float = 0.0,
        phase: float = 0.0,
        clip_shape: Part.Shape = None,
    ) -> Part.Shape:
        """
        Generate clipped slabs for one family at one phase.

        orientation_deg: in-plane rotation of the slab family (degrees)
        thickness:       slab width in the stacking direction (pre-clipping)
        tilt_deg:        out-of-plane tilt around the rib axis (degrees)
        phase:           fractional offset within one period (0 = rib, 0.5 = void)
        clip_shape:      shape to intersect with; defaults to self.shape
        """
        clip = clip_shape if clip_shape is not None else self.shape
        num_slabs = int(self.diag_xy / self.p.grid_spacing) + 3
        slabs = []

        for i in range(-num_slabs, num_slabs):
            slab = self._build_single_slab(
                index=i,
                phase=phase,
                thickness=thickness,
                orientation_deg=orientation_deg,
                tilt_deg=tilt_deg,
            )

            clipped = slab.common(clip)

            if clipped.isNull():
                continue

            if clipped.ShapeType == "Compound":
                solids = clipped.Solids
                if not solids:
                    continue
                clipped = Part.makeCompound(solids)

            slabs.append(clipped)

        if not slabs:
            return Part.makeCompound([])

        return Part.makeCompound(slabs)

    def _build_single_slab(
        self,
        index: int,
        phase: float,
        thickness: float,
        orientation_deg: float,
        tilt_deg: float,
    ) -> Part.Shape:
        """Build one unclipped slab box, positioned and rotated."""
        offset = (index + phase) * self.p.grid_spacing

        slab = Part.makeBox(
            thickness,
            self.ylen * 2,
            self.zlen * 2,
            FreeCAD.Vector(-thickness / 2, -self.ylen, -self.zlen),
        )

        slab.translate(FreeCAD.Vector(offset, 0, 0))
        slab.rotate(FreeCAD.Vector(0, 0, 0), FreeCAD.Vector(0, 0, 1), orientation_deg)

        if tilt_deg != 0.0:
            theta = m.radians(orientation_deg)
            tilt_axis = FreeCAD.Vector(-m.sin(theta), m.cos(theta), 0)
            slab.rotate(FreeCAD.Vector(0, 0, 0), tilt_axis, tilt_deg)

        slab.translate(self.center)
        return slab


# ============================================================
# FACE ANALYSIS
# ============================================================

class FaceAnalyzer:
    """Utility class for geometric face inspection."""

    @staticmethod
    def centroid(face) -> FreeCAD.Vector:
        return face.CenterOfMass

    @staticmethod
    def normal(face) -> FreeCAD.Vector:
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
      1 face  → extrude along face normal by `thickness`
      2 faces → shell between bottom (min Z centroid) and top (max Z centroid)

    Note: the 2-face shell path requires the two faces to share a compatible
    boundary for Part.makeShell to produce a closed solid. If the faces are
    non-planar or non-matching this will raise.
    """

    def __init__(self, shape: Part.Shape):
        self.shape = shape
        self.faces = shape.Faces
        if not self.faces:
            raise ValueError("Shape has no faces")

    def build(self, thickness: float = None) -> Part.Shape:
        if len(self.faces) == 1:
            return self._from_single_face(thickness)
        return self._from_two_faces()

    def _from_single_face(self, thickness: float) -> Part.Shape:
        face = self.faces[0]
        n = FaceAnalyzer.normal(face)
        t = thickness if thickness is not None else 1.0
        return face.extrude(n.multiply(t))

    def _from_two_faces(self) -> Part.Shape:
        sorted_faces = sorted(self.faces, key=lambda f: FaceAnalyzer.centroid(f).z)
        bottom = sorted_faces[0]
        top = sorted_faces[-1]
        # NOTE: makeShell from 2 open faces will only produce a valid closed
        # solid if the faces share a compatible boundary edge loop.
        shell = Part.makeShell([bottom, top])
        return Part.Solid(shell)


# ============================================================
# HIGH-LEVEL API WRAPPER
# ============================================================

class ThicknessGenerator:
    """User-facing interface for ThicknessBuilder."""

    def __init__(self, doc):
        self.doc = doc

    def from_object(
        self, obj_name: str, out_name: str = "ThickSolid", thickness: float = None
    ):
        obj = self.doc.getObject(obj_name)
        if obj is None:
            raise ValueError(f"Object '{obj_name}' not found in document")

        builder = ThicknessBuilder(obj.Shape)
        solid = builder.build(thickness=thickness)

        out = self.doc.addObject("Part::Feature", out_name)
        out.Shape = solid
        self.doc.recompute()
        return out


# ============================================================
# DIAGNOSTICS
# ============================================================

def see_objects(doc):
    for o in doc.Objects:
        print(o.Name, "|", o.Label, "|", o.TypeId)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    p = Params(
        grid_spacing=30,
        grid_angle=0,
        grid_orientation=0,
        grid_width=0.84,
        rib_clearance=2.0,   # 2mm gap between void slab and neighboring solid rib
    )

    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)

    see_objects(doc)

    full_wing = doc.getObject("Pad")
    bb = full_wing.Shape.BoundBox

    grid = TiltedGrid(bb=bb, params=p, shape=full_wing, doc=doc)
    grid.make_grid()

    doc.recompute()