import FreeCAD
import Part
import math

# ============================================================
# PARAMETERS
# ============================================================

# Maps plane name → (normal vector, two in-plane axes)
PLANE_DEFS = {
    'XY': {'normal': FreeCAD.Vector(0, 0, 1),
            'axis_u': FreeCAD.Vector(1, 0, 0),
            'axis_v': FreeCAD.Vector(0, 1, 0)},
    'XZ': {'normal': FreeCAD.Vector(0, 1, 0),
            'axis_u': FreeCAD.Vector(1, 0, 0),
            'axis_v': FreeCAD.Vector(0, 0, 1)},
    'YZ': {'normal': FreeCAD.Vector(1, 0, 0),
            'axis_u': FreeCAD.Vector(0, 1, 0),
            'axis_v': FreeCAD.Vector(0, 0, 1)},
}


class LWInfillParams:
    def __init__(self,
                 nozzle_diameter    = 0.4,
                 wall_thickness     = None,
                 rib_spacing        = 5.0,
                 rib_width          = None,
                 rib_angle          = 30.0,
                 grid_orientation   = 0.0,
                 primary_dir        = None,
                 construction_plane = 'XZ'):   # 'XY', 'XZ', or 'YZ'

        if construction_plane not in PLANE_DEFS:
            raise ValueError(f"construction_plane must be one of {list(PLANE_DEFS.keys())}")

        self.nozzle_diameter    = nozzle_diameter
        self.wall_thickness     = wall_thickness or nozzle_diameter
        self.rib_spacing        = rib_spacing
        self.rib_width          = rib_width or (2.0 * nozzle_diameter)
        self.rib_angle          = rib_angle
        self.grid_orientation   = grid_orientation
        self.construction_plane = construction_plane

        pdef = PLANE_DEFS[construction_plane]
        self.plane_normal = pdef['normal']
        self.plane_axis_u = pdef['axis_u']
        self.plane_axis_v = pdef['axis_v']

        # Project primary_dir onto the construction plane.
        # If result is null (e.g. Z passed for XY plane) fall back to auto-detect.
        self.primary_dir = self._project_primary(primary_dir)

    def _project_primary(self, pd):
        if pd is None:
            return None
        n = self.plane_normal
        dot = pd.x*n.x + pd.y*n.y + pd.z*n.z
        projected = FreeCAD.Vector(pd.x - dot*n.x,
                                   pd.y - dot*n.y,
                                   pd.z - dot*n.z)
        length = math.sqrt(projected.x**2 + projected.y**2 + projected.z**2)
        if length < 1e-6:
            print(f"[warn] primary_dir is parallel to {self.construction_plane} "
                  f"normal — using auto-detection.")
            return None
        return FreeCAD.Vector(projected.x/length,
                              projected.y/length,
                              projected.z/length)


# ============================================================
# HELPER: ROTATE VECTOR AROUND AN AXIS (Rodrigues)
# ============================================================
def rotate_vector_around(v, axis, angle_rad):
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    dot = v.x*axis.x + v.y*axis.y + v.z*axis.z
    cross = FreeCAD.Vector(
        axis.y*v.z - axis.z*v.y,
        axis.z*v.x - axis.x*v.z,
        axis.x*v.y - axis.y*v.x,
    )
    return FreeCAD.Vector(
        v.x*c + cross.x*s + axis.x*dot*(1-c),
        v.y*c + cross.y*s + axis.y*dot*(1-c),
        v.z*c + cross.z*s + axis.z*dot*(1-c),
    )


# ============================================================
# FACE SELECTION
# ============================================================
def find_face_on_plane(body, plane_normal):
    """
    Return the face whose centre of mass is most negative along
    the plane normal — i.e. the 'bottom' face for the chosen plane:
      XY (normal Z) → lowest Z face
      XZ (normal Y) → most negative Y face
      YZ (normal X) → most negative X face
    """
    best_face = None
    best_proj = float('inf')
    for face in body.Faces:
        try:
            com = face.CenterOfMass
            proj = com.x*plane_normal.x + com.y*plane_normal.y + com.z*plane_normal.z
            if proj < best_proj:
                best_proj = proj
                best_face = face
        except Exception:
            pass
    return best_face


# ============================================================
# GRID LINE GENERATION
# ============================================================
def create_angled_grid_lines(bb, params: LWInfillParams):
    """
    Generate two families of parallel lines in the construction plane.
    Works for XY, XZ, and YZ by using the plane axes from params.
    """
    n  = params.plane_normal
    au = params.plane_axis_u
    av = params.plane_axis_v

    corners = [
        FreeCAD.Vector(bb.XMin, bb.YMin, bb.ZMin),
        FreeCAD.Vector(bb.XMax, bb.YMin, bb.ZMin),
        FreeCAD.Vector(bb.XMin, bb.YMax, bb.ZMin),
        FreeCAD.Vector(bb.XMax, bb.YMax, bb.ZMin),
        FreeCAD.Vector(bb.XMin, bb.YMin, bb.ZMax),
        FreeCAD.Vector(bb.XMax, bb.YMin, bb.ZMax),
        FreeCAD.Vector(bb.XMin, bb.YMax, bb.ZMax),
        FreeCAD.Vector(bb.XMax, bb.YMax, bb.ZMax),
    ]

    def dot(v, a):
        return v.x*a.x + v.y*a.y + v.z*a.z

    span_u = max(dot(c, au) for c in corners) - min(dot(c, au) for c in corners)
    span_v = max(dot(c, av) for c in corners) - min(dot(c, av) for c in corners)

    # Primary direction within the construction plane
    if params.primary_dir is not None:
        pd = params.primary_dir
    else:
        pd = au if span_u >= span_v else av

    # In-plane perpendicular to primary = cross(normal, primary)
    perp_pd = FreeCAD.Vector(
        n.y*pd.z - n.z*pd.y,
        n.z*pd.x - n.x*pd.z,
        n.x*pd.y - n.y*pd.x,
    )
    perp_pd.normalize()

    ang = math.radians(params.rib_angle)
    rot = math.radians(params.grid_orientation)

    def make_family_dir(sign):
        d = FreeCAD.Vector(
            pd.x*math.cos(ang) + sign*perp_pd.x*math.sin(ang),
            pd.y*math.cos(ang) + sign*perp_pd.y*math.sin(ang),
            pd.z*math.cos(ang) + sign*perp_pd.z*math.sin(ang),
        )
        return rotate_vector_around(d, n, rot)

    d1 = make_family_dir(+1)
    d2 = make_family_dir(-1)

    center = FreeCAD.Vector(
        (bb.XMin + bb.XMax) / 2,
        (bb.YMin + bb.YMax) / 2,
        (bb.ZMin + bb.ZMax) / 2,
    )
    line_len = math.sqrt(bb.XLength**2 + bb.YLength**2 + bb.ZLength**2) * 2

    def generate_lines(d):
        # Stacking direction = in-plane perpendicular to d = cross(d, normal)
        stacking = FreeCAD.Vector(
            d.y*n.z - d.z*n.y,
            d.z*n.x - d.x*n.z,
            d.x*n.y - d.y*n.x,
        )
        stacking.normalize()

        proj_vals = [dot(c, stacking) for c in corners]
        min_p = min(proj_vals)
        max_p = max(proj_vals)
        num = int((max_p - min_p) / params.rib_spacing) + 3
        center_proj = dot(center, stacking)

        lines = []
        for i in range(-1, num + 1):
            offset = min_p + i * params.rib_spacing
            shift = offset - center_proj
            p0 = FreeCAD.Vector(
                center.x + stacking.x*shift,
                center.y + stacking.y*shift,
                center.z + stacking.z*shift,
            )
            start = FreeCAD.Vector(p0.x - d.x*line_len,
                                   p0.y - d.y*line_len,
                                   p0.z - d.z*line_len)
            end   = FreeCAD.Vector(p0.x + d.x*line_len,
                                   p0.y + d.y*line_len,
                                   p0.z + d.z*line_len)
            lines.append(Part.makeLine(start, end))
        return lines

    return generate_lines(d1), generate_lines(d2)


# ============================================================
# CREATE THIN RIB FACES FROM LINES
# ============================================================
def create_rib_faces(lines, plane_normal, rib_width):
    half_w = rib_width / 2.0
    faces = []
    for line in lines:
        try:
            start = line.Vertexes[0].Point
            end = line.Vertexes[-1].Point
            dir_vec = (end - start).normalize()
            perp = dir_vec.cross(plane_normal).normalize() * half_w
            p1 = start + perp
            p2 = start - perp
            p3 = end - perp
            p4 = end + perp
            wire = Part.makePolygon([p1, p2, p3, p4, p1])
            faces.append(Part.Face(wire))
        except Exception:
            pass
    return faces


# ============================================================
# MAIN GENERATION
# ============================================================
def generate_lw_infill(body, params: LWInfillParams, doc=None):
    if doc is None:
        doc = FreeCAD.ActiveDocument

    raw_shape = body.Shape if hasattr(body, "Shape") else body
    bb = raw_shape.BoundBox

    print(f"Construction plane: {params.construction_plane} "
          f"(normal: {params.plane_normal})")

    # 1. Find the reference face for this plane
    work_face = find_face_on_plane(raw_shape, params.plane_normal)
    if not work_face:
        raise ValueError("No valid face found on the body.")
    print(f"  Reference face area: {work_face.Area:.1f} mm²")

    # 2. Generate grid lines in the construction plane
    lines1, lines2 = create_angled_grid_lines(bb, params)
    print(f"  Family1 lines: {len(lines1)}, Family2 lines: {len(lines2)}")

    # 3. Convert to rib faces (thickness along plane_normal)
    rib_faces1 = create_rib_faces(lines1, params.plane_normal, params.rib_width)
    rib_faces2 = create_rib_faces(lines2, params.plane_normal, params.rib_width)
    print(f"  Rib faces: {len(rib_faces1)} + {len(rib_faces2)}")

    if not rib_faces1 and not rib_faces2:
        raise ValueError("No rib faces generated.")

    # 4. Extrude faces along plane_normal to span the full shape depth
    n = params.plane_normal

    def dot(v, a):
        return v.x*a.x + v.y*a.y + v.z*a.z

    corners = [
        FreeCAD.Vector(bb.XMin, bb.YMin, bb.ZMin),
        FreeCAD.Vector(bb.XMax, bb.YMax, bb.ZMax),
    ]
    min_along_n = min(dot(c, n) for c in corners)
    max_along_n = max(dot(c, n) for c in corners)
    extrude_depth = (max_along_n - min_along_n) * 1.2

    extrude_vec = FreeCAD.Vector(
        n.x * extrude_depth,
        n.y * extrude_depth,
        n.z * extrude_depth,
    )

    def extrude_faces(faces):
        solids = []
        for face in faces:
            try:
                face_proj = dot(face.CenterOfMass, n)
                shift = min_along_n - face_proj - extrude_depth * 0.1
                f = face.copy()
                f.translate(FreeCAD.Vector(n.x*shift, n.y*shift, n.z*shift))
                solids.append(f.extrude(extrude_vec))
            except Exception as e:
                print(f"  [warn] extrude failed: {e}")
        return solids

    solids1 = extrude_faces(rib_faces1)
    solids2 = extrude_faces(rib_faces2)
    print(f"  Extruded solids: {len(solids1)} + {len(solids2)}")

    def fuse_solids(solids):
        if not solids:
            return None
        fused = solids[0]
        for s in solids[1:]:
            try:
                fused = fused.fuse(s)
            except Exception:
                pass
        return fused

    fused1 = fuse_solids(solids1)
    fused2 = fuse_solids(solids2)

    rib_parts = [p for p in (fused1, fused2) if p is not None]
    if not rib_parts:
        raise ValueError("No rib solids to cut with.")
    fused_ribs = rib_parts[0]
    for p in rib_parts[1:]:
        try:
            fused_ribs = fused_ribs.fuse(p)
        except Exception:
            pass

    print("Cutting internal structure...")
    try:
        result = raw_shape.cut(fused_ribs)
    except Exception as e:
        print(f"Cut failed: {e}")
        raise

    obj = doc.addObject("Part::Feature", "WingWithInfill")
    obj.Shape = result
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

    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)
    see_objects(doc)

    wing = doc.getObject("Pad")
    if not wing:
        raise RuntimeError("Object 'Pad' not found in document.")

    params = LWInfillParams(
        nozzle_diameter  = 0.4,
        rib_spacing      = 10.0,
        rib_width        = 1.0,
        rib_angle        = 30.0,
        grid_orientation = 0.0,
        primary_dir      = FreeCAD.Vector(0, 0, 1),  # Z → auto-detected, no crash
        construction_plane = 'XZ'
    )

    result_body = generate_lw_infill(wing.Shape, params, doc)

    obj = doc.addObject("Part::Feature", "WingWithInfill")
    obj.Shape = result_body
    doc.recompute()
    print("Lightweight wing created successfully.")