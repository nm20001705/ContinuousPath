import FreeCAD
import Part
import math as m


# ============================================================
# PARAMETERS
# ============================================================

class LWInfillParams:
    """
    Parameters for a two-family tilted slab lattice with elliptical
    lightening holes, designed for LW-PLA vase-mode printing.

    nozzle_diameter:
        Nozzle width (mm). Used to derive rib_width if not set explicitly.

    rib_spacing:
        Center-to-center distance between ribs within one family (mm).

    rib_width:
        Thickness of each rib slab (mm). Defaults to 2 * nozzle_diameter.

    rib_angle:
        Half-angle between the two rib families (degrees).
        Family 1 runs at +rib_angle from primary_dir,
        family 2 runs at -rib_angle from primary_dir.

    tilt_deg:
        Out-of-plane tilt of each family around its own rib axis (degrees).
        Family 1 tilts +tilt_deg, family 2 tilts -tilt_deg.
        This makes the crossing vector horizontal (in XY) rather than vertical.

    grid_orientation:
        Global rotation of both families around Z (degrees).

    primary_dir:
        Optional span direction (FreeCAD.Vector in XY).
        If None, the longer bounding box axis is used.

    hole_margin:
        Minimum material between ellipse edge and rib boundary (mm).
        Enforced in the canonical frame before tilting.

    min_hole_size:
        Skip lightening if both semi-axes are smaller than this (mm).
    """

    def __init__(
        self,
        nozzle_diameter  = 0.4,
        rib_spacing      = 10.0,
        rib_width        = None,
        rib_angle        = 30.0,
        tilt_deg         = 30.0,
        grid_orientation = 0.0,
        primary_dir      = None,
        hole_margin      = 1.0,
        min_hole_size    = 4.0,
    ):
        self.nozzle_diameter  = nozzle_diameter
        self.rib_spacing      = rib_spacing
        self.rib_width        = rib_width or (2.0 * nozzle_diameter)
        self.rib_angle        = rib_angle
        self.tilt_deg         = tilt_deg
        self.grid_orientation = grid_orientation
        self.primary_dir      = primary_dir
        self.hole_margin      = hole_margin
        self.min_hole_size    = min_hole_size

    @property
    def semi_major(self) -> float:
        """
        Half-length of the ellipse along the rib (X in canonical frame).
        Capped at rib_spacing/2 - hole_margin so holes never reach node regions.
        """
        return self.rib_spacing / 2.0 - self.hole_margin

    @property
    def semi_minor(self) -> float:
        """
        Half-width of the ellipse across the rib (Y in canonical frame).
        Capped at rib_width/2 - hole_margin so holes never reach slab edges.
        """
        return self.rib_width / 2.0 - self.hole_margin

    def validate(self):
        if self.semi_minor <= 0:
            raise ValueError(
                f"hole_margin ({self.hole_margin}) too large for "
                f"rib_width ({self.rib_width}): semi_minor = {self.semi_minor:.3f}."
            )
        if self.semi_major <= 0:
            raise ValueError(
                f"hole_margin ({self.hole_margin}) too large for "
                f"rib_spacing ({self.rib_spacing}): semi_major = {self.semi_major:.3f}."
            )


# ============================================================
# HELPERS
# ============================================================

def is_valid(shape) -> bool:
    if shape is None or shape.isNull():
        return False
    try:
        return len(shape.Solids) > 0
    except Exception:
        return False


def safe_fuse(a: Part.Shape, b: Part.Shape) -> Part.Shape:
    if not is_valid(a):
        return b
    if not is_valid(b):
        return a
    try:
        return a.fuse(b)
    except Exception as e:
        print(f"[warn] fuse failed, using compound: {e}")
        return Part.makeCompound(a.Solids + b.Solids)


# ============================================================
# FAMILY DIRECTION COMPUTATION
# ============================================================

def compute_family_directions(bb, params: LWInfillParams):
    """
    Compute the rib direction vectors for both families.

    Returns (dir1, dir2) — unit vectors in XY for each family.

    The primary direction is either user-supplied or derived from the
    longer bounding box axis. Each family is rotated ±rib_angle from
    primary, then the whole grid is rotated by grid_orientation.
    """
    # Primary direction
    if params.primary_dir is not None:
        pd = FreeCAD.Vector(params.primary_dir.x, params.primary_dir.y, 0)
        pd.normalize()
    else:
        if bb.XLength >= bb.YLength:
            pd = FreeCAD.Vector(1, 0, 0)
        else:
            pd = FreeCAD.Vector(0, 1, 0)

    # Perpendicular to primary in XY
    perp = FreeCAD.Vector(-pd.y, pd.x, 0)

    ang = m.radians(params.rib_angle)
    rot = m.radians(params.grid_orientation)

    def rotate_xy(v, angle):
        return FreeCAD.Vector(
            v.x * m.cos(angle) - v.y * m.sin(angle),
            v.x * m.sin(angle) + v.y * m.cos(angle),
            0,
        )

    # Family directions before global rotation
    d1_local = FreeCAD.Vector(
        pd.x * m.cos(ang) + perp.x * m.sin(ang),
        pd.y * m.cos(ang) + perp.y * m.sin(ang),
        0,
    )
    d2_local = FreeCAD.Vector(
        pd.x * m.cos(ang) - perp.x * m.sin(ang),
        pd.y * m.cos(ang) - perp.y * m.sin(ang),
        0,
    )

    dir1 = rotate_xy(d1_local, rot)
    dir2 = rotate_xy(d2_local, rot)
    dir1.normalize()
    dir2.normalize()

    return dir1, dir2


# ============================================================
# SLAB BUILDER
# ============================================================

class SlabFamily:
    """
    Builds one family of tilted lightened slabs.

    In the canonical frame:
      - The slab runs along X (the rib direction)
      - The slab cross-section is rib_width in Y, oversized in Z
      - Elliptical holes are cut periodically along X before any rotation

    Then the slab is:
      1. Rotated in XY to align with the family direction
      2. Tilted around its own rib axis by tilt_deg
      3. Translated to the correct position
      4. Clipped to the wing shape
    """

    def __init__(self, shape: Part.Shape, params: LWInfillParams,
                 direction: FreeCAD.Vector, tilt_deg: float):
        self.shape     = shape
        self.p         = params
        self.direction = direction   # unit vector in XY
        self.tilt_deg  = tilt_deg

        bb = shape.BoundBox
        self.bb      = bb
        self.xlen    = bb.XLength
        self.ylen    = bb.YLength
        self.zlen    = bb.ZLength
        self.diag_xy = m.sqrt(self.xlen**2 + self.ylen**2)
        self.diag_3d = m.sqrt(self.xlen**2 + self.ylen**2 + self.zlen**2)

        self.center = FreeCAD.Vector(
            bb.XMin + self.xlen / 2,
            bb.YMin + self.ylen / 2,
            bb.ZMin + self.zlen / 2,
        )

        # Orientation angle of this family (degrees, measured from X axis)
        self.orientation_deg = m.degrees(m.atan2(direction.y, direction.x))

        # Tilt axis = in-plane perpendicular to the rib direction
        self.tilt_axis = FreeCAD.Vector(-direction.y, direction.x, 0)

    def build(self) -> Part.Shape:
        """Build all slabs for this family, clipped to the wing shape."""
        # Project bounding box corners onto the family direction to get
        # tight iteration bounds (avoids generating slabs far outside shape)
        corners = [
            FreeCAD.Vector(self.bb.XMin, self.bb.YMin, 0),
            FreeCAD.Vector(self.bb.XMax, self.bb.YMin, 0),
            FreeCAD.Vector(self.bb.XMin, self.bb.YMax, 0),
            FreeCAD.Vector(self.bb.XMax, self.bb.YMax, 0),
        ]
        perp = self.tilt_axis  # perpendicular to rib in XY = stacking direction
        proj = [c.dot(perp) for c in corners]
        center_proj = self.center.dot(perp)
        span = max(proj) - min(proj)
        num = int(span / self.p.rib_spacing) + 3

        slabs = []
        for i in range(-num, num + 1):
            slab = self._build_one_slab(i)
            clipped = slab.common(self.shape)
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

    def _build_one_slab(self, index: int) -> Part.Shape:
        """
        Build one lightened slab in canonical frame, then orient it.

        Canonical frame: slab along X, cross-section in Y*Z, centered at origin.
        Stacking offset is applied along the perpendicular (tilt_axis) direction.
        """
        w      = self.p.rib_width
        length = self.diag_3d * 2   # oversized — clipped later

        # ---- 1. Canonical slab box ----
        slab = Part.makeBox(
            length, w, self.zlen * 2,
            FreeCAD.Vector(-length / 2, -w / 2, -self.zlen),
        )

        # ---- 2. Cut periodic elliptical holes in canonical frame ----
        if (self.p.semi_minor >= self.p.min_hole_size or
                self.p.semi_major >= self.p.min_hole_size):
            holes = self._build_hole_strip(length)
            if is_valid(holes):
                slab = slab.cut(holes)

        # ---- 3. Rotate in XY to align with family direction ----
        slab.rotate(
            FreeCAD.Vector(0, 0, 0),
            FreeCAD.Vector(0, 0, 1),
            self.orientation_deg,
        )

        # ---- 4. Tilt around the rib axis ----
        if self.tilt_deg != 0.0:
            slab.rotate(
                FreeCAD.Vector(0, 0, 0),
                self.tilt_axis,
                self.tilt_deg,
            )

        # ---- 5. Translate: stack offset along tilt_axis + shape center ----
        offset_vec = FreeCAD.Vector(
            self.tilt_axis.x * index * self.p.rib_spacing + self.center.x,
            self.tilt_axis.y * index * self.p.rib_spacing + self.center.y,
            self.center.z,
        )
        slab.translate(offset_vec)

        return slab

    def _build_hole_strip(self, length: float) -> Part.Shape:
        """
        Periodic elliptical cylinders along X in the canonical frame.
        One hole per rib_spacing period, extruded in Z to fully penetrate.
        """
        rx = self.p.semi_major
        ry = self.p.semi_minor
        s  = self.p.rib_spacing
        z_depth = self.zlen * 2 + 2.0

        num_holes = int(length / s) + 2
        holes = []

        for k in range(-num_holes, num_holes + 1):
            cx = k * s
            ellipse = Part.Ellipse(FreeCAD.Vector(0, 0, 0), rx, ry)
            wire = Part.Wire(ellipse.toShape())
            face = Part.Face(wire)
            face.translate(FreeCAD.Vector(cx, 0, -z_depth / 2))
            cyl = face.extrude(FreeCAD.Vector(0, 0, z_depth))
            holes.append(cyl)

        if not holes:
            return Part.makeCompound([])
        return Part.makeCompound(holes)


# ============================================================
# MAIN GENERATOR
# ============================================================

def generate_lw_infill(shape, params: LWInfillParams, doc=None):
    """
    Generate a two-family tilted slab lattice with elliptical lightening,
    clipped to `shape`, and add the result to `doc`.
    """
    if doc is None:
        doc = FreeCAD.ActiveDocument

    params.validate()

    raw_shape = shape.Shape if hasattr(shape, "Shape") else shape
    bb = raw_shape.BoundBox

    # Compute family directions
    dir1, dir2 = compute_family_directions(bb, params)
    print(f"Family 1 direction: ({dir1.x:.3f}, {dir1.y:.3f})")
    print(f"Family 2 direction: ({dir2.x:.3f}, {dir2.y:.3f})")

    # Build each family
    print("Generating family 1...")
    fam1 = SlabFamily(raw_shape, params, dir1, +params.tilt_deg)
    ribs1 = fam1.build()
    print(f"  {len(ribs1.Solids) if is_valid(ribs1) else 0} slab segments")

    print("Generating family 2...")
    fam2 = SlabFamily(raw_shape, params, dir2, -params.tilt_deg)
    ribs2 = fam2.build()
    print(f"  {len(ribs2.Solids) if is_valid(ribs2) else 0} slab segments")

    # Combine — use makeCompound to avoid expensive boolean fuse
    parts = [s for s in [ribs1, ribs2] if is_valid(s)]
    if not parts:
        raise ValueError("No rib geometry produced.")

    final = Part.makeCompound(parts)

    obj = doc.addObject("Part::Feature", "LWInfill")
    obj.Shape = final
    doc.recompute()
    print("Done.")
    return obj


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

    params = LWInfillParams(
        nozzle_diameter  = 0.4,
        rib_spacing      = 20.0,
        rib_width        = 3.0,
        rib_angle        = 30.0,    # ±30° from primary direction
        tilt_deg         = 30.0,    # tilt makes crossing vector horizontal
        grid_orientation = 0.0,
        primary_dir      = None,    # auto-detect from bounding box
        hole_margin      = 0.5,
        min_hole_size    = 4.0,
    )

    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)
    see_objects(doc)

    wing = doc.getObject("Pad")
    generate_lw_infill(wing, params, doc)

    doc.recompute()