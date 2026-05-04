import FreeCAD
import Part
import math

# ============================================================
# PARAMETERS
# ============================================================

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
                 construction_plane = 'XZ'):
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
            print(f"[warn] primary_dir is parallel to {self.construction_plane} normal — using auto-detection.")
            return None
        return FreeCAD.Vector(projected.x/length, projected.y/length, projected.z/length)


def rotate_vector_around(v, axis, angle_rad):
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    dot_v = v.x*axis.x + v.y*axis.y + v.z*axis.z
    cross = FreeCAD.Vector(axis.y*v.z - axis.z*v.y,
                           axis.z*v.x - axis.x*v.z,
                           axis.x*v.y - axis.y*v.x)
    return FreeCAD.Vector(
        v.x*c + cross.x*s + axis.x*dot_v*(1-c),
        v.y*c + cross.y*s + axis.y*dot_v*(1-c),
        v.z*c + cross.z*s + axis.z*dot_v*(1-c)
    )


def find_face_on_plane(body, plane_normal):
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


def create_angled_grid_lines(bb, params):
    n = params.plane_normal
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

    if params.primary_dir is not None:
        pd = params.primary_dir
    else:
        pd = au if span_u >= span_v else av

    perp_pd = FreeCAD.Vector(n.y*pd.z - n.z*pd.y,
                             n.z*pd.x - n.x*pd.z,
                             n.x*pd.y - n.y*pd.x)
    perp_pd.normalize()

    ang = math.radians(params.rib_angle)
    rot = math.radians(params.grid_orientation)

    def make_family_dir(sign):
        d = FreeCAD.Vector(pd.x*math.cos(ang) + sign*perp_pd.x*math.sin(ang),
                           pd.y*math.cos(ang) + sign*perp_pd.y*math.sin(ang),
                           pd.z*math.cos(ang) + sign*perp_pd.z*math.sin(ang))
        return rotate_vector_around(d, n, rot)

    d1 = make_family_dir(+1)
    d2 = make_family_dir(-1)

    center = FreeCAD.Vector((bb.XMin+bb.XMax)/2,
                            (bb.YMin+bb.YMax)/2,
                            (bb.ZMin+bb.ZMax)/2)
    line_len = math.sqrt(bb.XLength**2 + bb.YLength**2 + bb.ZLength**2) * 2

    def generate_lines(d):
        stacking = FreeCAD.Vector(d.y*n.z - d.z*n.y,
                                  d.z*n.x - d.x*n.z,
                                  d.x*n.y - d.y*n.x)
        stacking.normalize()
        proj_vals = [dot(c, stacking) for c in corners]
        min_p = min(proj_vals)
        max_p = max(proj_vals)
        num = int((max_p - min_p) / params.rib_spacing) + 3
        center_proj = dot(center, stacking)

        lines = []
        for i in range(-1, num+1):
            offset = min_p + i * params.rib_spacing
            shift = offset - center_proj
            p0 = FreeCAD.Vector(center.x + stacking.x*shift,
                                center.y + stacking.y*shift,
                                center.z + stacking.z*shift)
            start = FreeCAD.Vector(p0.x - d.x*line_len,
                                   p0.y - d.y*line_len,
                                   p0.z - d.z*line_len)
            end   = FreeCAD.Vector(p0.x + d.x*line_len,
                                   p0.y + d.y*line_len,
                                   p0.z + d.z*line_len)
            lines.append(Part.makeLine(start, end))
        return lines

    return generate_lines(d1), generate_lines(d2)


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


def add_bridges_along_ribs(cut_body, rib_center_lines, rib_width, plane_normal, bridge_size=None):
    if bridge_size is None:
        bridge_size = rib_width * 0.5
    if cut_body.ShapeType in ("Compound", "CompSolid"):
        solids = cut_body.Solids
    else:
        solids = [cut_body]

    n = plane_normal.normalize()
    bridges = []

    for line in rib_center_lines:
        start = line.Vertexes[0].Point
        end = line.Vertexes[-1].Point
        mid = (start + end) * 0.5
        rib_dir = (end - start).normalize()
        span_dir = n
        third_dir = rib_dir.cross(span_dir).normalize()

        len_along_rib = bridge_size
        len_along_span = rib_width * 1.1
        len_along_third = bridge_size

        box = Part.makeBox(len_along_rib, len_along_span, len_along_third)
        box.translate(FreeCAD.Vector(-len_along_rib/2, -len_along_span/2, -len_along_third/2))

        mat = FreeCAD.Matrix()
        mat.A11 = rib_dir.x;   mat.A12 = span_dir.x;   mat.A13 = third_dir.x
        mat.A21 = rib_dir.y;   mat.A22 = span_dir.y;   mat.A23 = third_dir.y
        mat.A31 = rib_dir.z;   mat.A32 = span_dir.z;   mat.A33 = third_dir.z
        mat.A44 = 1.0
        rot = FreeCAD.Rotation(mat)
        box.Placement = FreeCAD.Placement(mid, rot)
        bridges.append(box)

    if not bridges:
        return cut_body
    all_parts = solids + bridges
    unified = all_parts[0]
    for p in all_parts[1:]:
        try:
            unified = unified.fuse(p)
        except Exception as e:
            print(f"  [warn] fuse failed: {e}")
    return unified


def create_cut_result(body, params, doc=None):
    if doc is None:
        doc = FreeCAD.ActiveDocument
    raw_shape = body.Shape if hasattr(body, "Shape") else body
    bb = raw_shape.BoundBox

    lines1, lines2 = create_angled_grid_lines(bb, params)
    rib_faces1 = create_rib_faces(lines1, params.plane_normal, params.rib_width)
    rib_faces2 = create_rib_faces(lines2, params.plane_normal, params.rib_width)

    n = params.plane_normal
    def dot(v, a):
        return v.x*a.x + v.y*a.y + v.z*a.z

    corners = [FreeCAD.Vector(bb.XMin, bb.YMin, bb.ZMin),
               FreeCAD.Vector(bb.XMax, bb.YMax, bb.ZMax)]
    min_along_n = min(dot(c, n) for c in corners)
    max_along_n = max(dot(c, n) for c in corners)
    extrude_depth = (max_along_n - min_along_n) * 1.2
    extrude_vec = FreeCAD.Vector(n.x*extrude_depth, n.y*extrude_depth, n.z*extrude_depth)

    def extrude_faces(faces):
        solids = []
        for face in faces:
            try:
                face_proj = dot(face.CenterOfMass, n)
                shift = min_along_n - face_proj - extrude_depth * 0.1
                f = face.copy()
                f.translate(FreeCAD.Vector(n.x*shift, n.y*shift, n.z*shift))
                solids.append(f.extrude(extrude_vec))
            except Exception:
                pass
        return solids

    solids1 = extrude_faces(rib_faces1)
    solids2 = extrude_faces(rib_faces2)

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
        raise ValueError("No rib solids generated.")
    fused_ribs = rib_parts[0]
    for p in rib_parts[1:]:
        try:
            fused_ribs = fused_ribs.fuse(p)
        except Exception:
            pass

    cut_result = raw_shape.cut(fused_ribs)
    all_center_lines = lines1 + lines2
    return cut_result, all_center_lines


def generate_lw_infill(body, params, doc=None):
    cut_result, center_lines = create_cut_result(body, params, doc)
    bridged_result = add_bridges_along_ribs(cut_result, center_lines,
                                            params.rib_width, params.plane_normal,
                                            bridge_size=params.rib_width * 10)
    obj = doc.addObject("Part::Feature", "WingWithInfill")
    obj.Shape = bridged_result
    doc.recompute()
    return obj


# ============================================================
# NEW: POINT COLLECTION AND RIB ASSIGNMENT
# ============================================================
def collect_midpoints_per_rib(cut_body, rib_center_lines, z_min, z_max, z_step, max_dist=5.0):
    """
    Slice cut_body horizontally, collect midpoints of intersection edges,
    assign each midpoint to the nearest rib centre line.
    """
    ribs = list(enumerate(rib_center_lines))
    points_by_rib = {idx: [] for idx, _ in ribs}

    # Get list of solids (handle compound or single shape)
    if hasattr(cut_body, 'Solids') and cut_body.Solids:
        solids = cut_body.Solids
    else:
        solids = [cut_body]

    def process_edge(edge, points_dict):
        if len(edge.Vertexes) >= 2:
            p1 = edge.Vertexes[0].Point
            p2 = edge.Vertexes[-1].Point
            mid = (p1 + p2) * 0.5
            best_idx = None
            best_dist = float('inf')
            for idx, line in ribs:
                start = line.Vertexes[0].Point
                end = line.Vertexes[-1].Point
                line_dir = (end - start).normalize()
                t = (mid - start).dot(line_dir)
                if t < 0:
                    dist = mid.distanceToPoint(start)
                elif t > (end - start).Length:
                    dist = mid.distanceToPoint(end)
                else:
                    proj = start + line_dir * t
                    dist = mid.distanceToPoint(proj)
                if dist < best_dist:
                    best_dist = dist
                    best_idx = idx
            if best_idx is not None and best_dist <= max_dist:
                points_dict[best_idx].append(mid)

    z = z_min
    while z <= z_max:
        direction = FreeCAD.Vector(0, 0, 1)
        distance = z
        for solid in solids:
            try:
                slice_result = solid.slice(direction, distance)
                # slice_result can be a Compound or a list of Compounds
                if hasattr(slice_result, 'Edges'):
                    for edge in slice_result.Edges:
                        process_edge(edge, points_by_rib)
                elif isinstance(slice_result, list):
                    for item in slice_result:
                        if hasattr(item, 'Edges'):
                            for edge in item.Edges:
                                process_edge(edge, points_by_rib)
                else:
                    # fallback: try to iterate edges directly
                    try:
                        for edge in slice_result.Edges:
                            process_edge(edge, points_by_rib)
                    except:
                        pass
            except Exception:
                # fallback to section
                try:
                    section = solid.section(FreeCAD.Vector(0,0,z), direction)
                    for edge in section.Edges:
                        process_edge(edge, points_by_rib)
                except Exception:
                    continue
        z += z_step
    return points_by_rib


def show_points_per_rib(points_by_rib, doc):
    """Create one point cloud object per rib."""
    count = 0
    for idx, pts in points_by_rib.items():
        if not pts:
            continue
        vertices = [Part.Vertex(p) for p in pts]
        compound = Part.Compound(vertices)
        obj = doc.addObject("Part::Feature", f"RibPoints_{idx}")
        obj.Shape = compound
        count += 1
    doc.recompute()
    print(f"Visualized {count} ribs with points (total points = {sum(len(p) for p in points_by_rib.values())})")


# ============================================================
# MAIN (single, clean block)
# ============================================================
if __name__ == "__main__":
    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)

    wing = doc.getObject("Pad")
    if not wing:
        raise RuntimeError("Object 'Pad' not found in document.")

    # Use smaller rib_spacing to ensure ribs exist
    params = LWInfillParams(
        nozzle_diameter  = 0.4,
        rib_spacing      = 70.0,          # smaller than 40 to get ribs
        rib_width        = 0.1,
        rib_angle        = 30.0,
        grid_orientation = 0.0,
        primary_dir      = FreeCAD.Vector(0, 0, 1),
        construction_plane = 'XZ'
    )

    # 1. Create cut result (wing minus ribs)
    cut_result, all_center_lines = create_cut_result(wing.Shape, params, doc)
    print("Cut created.")

    # 2. Visualise original rib centre lines
    rib_lines_compound = Part.Compound([line for line in all_center_lines])
    obj_lines = doc.addObject("Part::Feature", "RibCentreLines")
    obj_lines.Shape = rib_lines_compound
    doc.recompute()
    print(f"Visualized {len(all_center_lines)} rib centre lines.")

    # 3. Slice range
    bb = cut_result.BoundBox
    if bb.ZMax - bb.ZMin < 1.0:
        print("Wing height too small – cannot slice.")
    else:
        z_min = bb.ZMin + 0.5
        z_max = bb.ZMax - 0.5
        z_step = 10.0
        print(f"Slicing from z={z_min:.1f} to {z_max:.1f} step {z_step:.1f}")

        # 4. Collect and assign midpoints
        points_by_rib = collect_midpoints_per_rib(cut_result, all_center_lines,
                                                  z_min, z_max, z_step, max_dist=5.0)

        # 5. Show point clouds per rib
        show_points_per_rib(points_by_rib, doc)

    # (Optional) Generate final solid with bridges
    # result_body = generate_lw_infill(wing.Shape, params, doc)
    # print("Lightweight wing with bridges created successfully.")