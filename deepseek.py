import FreeCAD
import Part
import math
import itertools

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

        # Build rotation matrix
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
# MIDPOINT COLLECTION AND CURVE VISUALIZATION
# ============================================================
def collect_midpoints_from_slices(cut_body, z_min, z_max, z_step):
    if cut_body.ShapeType in ("Compound", "CompSolid"):
        solids = cut_body.Solids
    else:
        solids = [cut_body]

    all_points = []
    z = z_min
    while z <= z_max:
        plane_origin = FreeCAD.Vector(0, 0, z)
        plane_normal = FreeCAD.Vector(0, 0, 1)
        for solid in solids:
            try:
                section = solid.section(plane_origin, plane_normal)
                for edge in section.Edges:
                    if len(edge.Vertexes) >= 2:
                        p1 = edge.Vertexes[0].Point
                        p2 = edge.Vertexes[-1].Point
                        all_points.append((p1 + p2) * 0.5)
            except Exception:
                continue
        z += z_step
    return all_points


def cluster_points_by_proximity(points, distance_threshold=5.0):
    """
    Simple clustering: start with first unassigned point, collect all points within distance.
    Returns list of clusters (each a list of points).
    """
    if not points:
        return []
    points_copy = list(points)
    clusters = []
    while points_copy:
        seed = points_copy.pop(0)
        cluster = [seed]
        i = 0
        while i < len(points_copy):
            p = points_copy[i]
            if min(p.distanceToPoint(q) for q in cluster) < distance_threshold:
                cluster.append(points_copy.pop(i))
            else:
                i += 1
        clusters.append(cluster)
    return clusters


def create_curves_from_clusters(clusters, doc):
    """Create a separate compound of wires/curves for all clusters."""
    curves = []
    for cluster in clusters:
        if len(cluster) < 2:
            continue
        # Sort by Z coordinate (since slices are horizontal)
        cluster_sorted = sorted(cluster, key=lambda p: p.z)
        # Create a polyline (or use BSpline for smoothness)
        polyline = Part.makePolygon(cluster_sorted)
        # The last point repeats the first, so we need to drop the extra
        if len(cluster_sorted) > 2:
            # Remove the final closing edge
            edges = polyline.Edges[:-1]
            wire = Part.Wire(edges)
        else:
            wire = Part.Wire(polyline)
        curves.append(wire)

    if not curves:
        return None
    compound = Part.Compound(curves)
    obj = doc.addObject("Part::Feature", "RibCentreLines")
    obj.Shape = compound
    doc.recompute()
    print(f"Created {len(curves)} rib centre curves (total {sum(len(c) for c in clusters)} points)")
    return obj


def visualize_rib_centre_lines(cut_body, z_step=0.5, distance_threshold=5.0, doc=None):
    bb = cut_body.BoundBox
    z_min = bb.ZMin + 0.1
    z_max = bb.ZMax - 0.1
    points = collect_midpoints_from_slices(cut_body, z_min, z_max, z_step)
    if not points:
        print("No midpoints collected.")
        return
    clusters = cluster_points_by_proximity(points, distance_threshold)
    create_curves_from_clusters(clusters, doc)


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)

    wing = doc.getObject("Pad")
    if not wing:
        raise RuntimeError("Object 'Pad' not found in document.")

    params = LWInfillParams(
        nozzle_diameter  = 0.4,
        rib_spacing      = 40.0,
        rib_width        = 0.1,
        rib_angle        = 30.0,
        grid_orientation = 0.0,
        primary_dir      = FreeCAD.Vector(0, 0, 1),
        construction_plane = 'XZ'
    )

    # --- Get rib centre lines (directly from generation) ---
    bb = wing.Shape.BoundBox
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_lines = lines1 + lines2

    # --- Create a compound of these lines and add to document as a separate object ---
    compounds = [line for line in all_lines]  # list of edges
    rib_lines_compound = Part.Compound(compounds)
    obj_lines = doc.addObject("Part::Feature", "RibCentreLines")
    obj_lines.Shape = rib_lines_compound
    doc.recompute()
    print(f"Visualized {len(all_lines)} rib centre lines.")

    # --- Now also generate the final solid with bridges if desired ---
    result_body = generate_lw_infill(wing.Shape, params, doc)
    print("Lightweight wing with bridges created successfully.")