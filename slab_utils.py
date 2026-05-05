# slab_utils.py
import FreeCAD
import Part
import math

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
    def __init__(self, nozzle_diameter=0.4, wall_thickness=None,
                 rib_spacing=5.0, rib_width=None, rib_angle=30.0,
                 grid_orientation=0.0, primary_dir=None,
                 construction_plane='XZ'):
        if construction_plane not in PLANE_DEFS:
            raise ValueError(f"construction_plane must be one of {list(PLANE_DEFS.keys())}")
        self.nozzle_diameter = nozzle_diameter
        self.wall_thickness = wall_thickness or nozzle_diameter
        self.rib_spacing = rib_spacing
        self.rib_width = rib_width or (2.0 * nozzle_diameter)
        self.rib_angle = rib_angle
        self.grid_orientation = grid_orientation
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
        proj = FreeCAD.Vector(pd.x - dot*n.x, pd.y - dot*n.y, pd.z - dot*n.z)
        L = proj.Length
        if L < 1e-6:
            print(f"[warn] primary_dir parallel to {self.construction_plane} normal – auto-detect.")
            return None
        return proj.normalize()

def rotate_vector_around(v, axis, angle_rad):
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    dot = v.x*axis.x + v.y*axis.y + v.z*axis.z
    cross = FreeCAD.Vector(axis.y*v.z - axis.z*v.y,
                           axis.z*v.x - axis.x*v.z,
                           axis.x*v.y - axis.y*v.x)
    return FreeCAD.Vector(
        v.x*c + cross.x*s + axis.x*dot*(1-c),
        v.y*c + cross.y*s + axis.y*dot*(1-c),
        v.z*c + cross.z*s + axis.z*dot*(1-c)
    )

def create_angled_grid_lines(bb, params):
    n, au, av = params.plane_normal, params.plane_axis_u, params.plane_axis_v
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
    def dot(v, a): return v.x*a.x + v.y*a.y + v.z*a.z
    span_u = max(dot(c, au) for c in corners) - min(dot(c, au) for c in corners)
    span_v = max(dot(c, av) for c in corners) - min(dot(c, av) for c in corners)
    if params.primary_dir:
        pd = params.primary_dir
    else:
        pd = au if span_u >= span_v else av
    perp_pd = FreeCAD.Vector(n.y*pd.z - n.z*pd.y,
                             n.z*pd.x - n.x*pd.z,
                             n.x*pd.y - n.y*pd.x).normalize()
    ang = math.radians(params.rib_angle)
    rot = math.radians(params.grid_orientation)
    def make_family_dir(sign):
        d = FreeCAD.Vector(pd.x*math.cos(ang) + sign*perp_pd.x*math.sin(ang),
                           pd.y*math.cos(ang) + sign*perp_pd.y*math.sin(ang),
                           pd.z*math.cos(ang) + sign*perp_pd.z*math.sin(ang))
        return rotate_vector_around(d, n, rot)
    d1, d2 = make_family_dir(1), make_family_dir(-1)
    center = FreeCAD.Vector((bb.XMin+bb.XMax)/2, (bb.YMin+bb.YMax)/2, (bb.ZMin+bb.ZMax)/2)
    line_len = math.sqrt(bb.XLength**2 + bb.YLength**2 + bb.ZLength**2) * 2
    def generate_lines(d):
        stacking = FreeCAD.Vector(d.y*n.z - d.z*n.y, d.z*n.x - d.x*n.z, d.x*n.y - d.y*n.x).normalize()
        proj_vals = [dot(c, stacking) for c in corners]
        min_p, max_p = min(proj_vals), max(proj_vals)
        num = int((max_p - min_p) / params.rib_spacing) + 3
        center_proj = dot(center, stacking)
        lines = []
        for i in range(-1, num+1):
            offset = min_p + i * params.rib_spacing
            shift = offset - center_proj
            p0 = FreeCAD.Vector(center.x + stacking.x*shift,
                                center.y + stacking.y*shift,
                                center.z + stacking.z*shift)
            start = FreeCAD.Vector(p0.x - d.x*line_len, p0.y - d.y*line_len, p0.z - d.z*line_len)
            end   = FreeCAD.Vector(p0.x + d.x*line_len, p0.y + d.y*line_len, p0.z + d.z*line_len)
            lines.append(Part.makeLine(start, end))
        return lines
    return generate_lines(d1), generate_lines(d2)

def create_rib_faces(lines, plane_normal, rib_width):
    half_w = rib_width / 2.0
    faces = []
    for line in lines:
        try:
            start, end = line.Vertexes[0].Point, line.Vertexes[-1].Point
            dir_vec = (end - start).normalize()
            perp = dir_vec.cross(plane_normal).normalize() * half_w
            p1 = start + perp
            p2 = start - perp
            p3 = end - perp
            p4 = end + perp
            wire = Part.makePolygon([p1, p2, p3, p4, p1])
            faces.append(Part.Face(wire))
        except:
            pass
    return faces

def extrude_rib_faces_to_solids(faces, plane_normal, bb):
    n = plane_normal.normalize()
    def dot(v, a): return v.x*a.x + v.y*a.y + v.z*a.z
    corners = [FreeCAD.Vector(bb.XMin, bb.YMin, bb.ZMin),
               FreeCAD.Vector(bb.XMax, bb.YMax, bb.ZMax)]
    min_along_n = min(dot(c, n) for c in corners)
    max_along_n = max(dot(c, n) for c in corners)
    extrude_depth = (max_along_n - min_along_n) * 1.2
    extrude_vec = n * extrude_depth
    solids = []
    for face in faces:
        try:
            face_proj = dot(face.CenterOfMass, n)
            shift = min_along_n - face_proj - extrude_depth * 0.1
            f = face.copy()
            f.translate(n * shift)
            solids.append(f.extrude(extrude_vec))
        except:
            continue
    return solids

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
        start, end = line.Vertexes[0].Point, line.Vertexes[-1].Point
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
        mat.A11, mat.A12, mat.A13 = rib_dir.x, span_dir.x, third_dir.x
        mat.A21, mat.A22, mat.A23 = rib_dir.y, span_dir.y, third_dir.y
        mat.A31, mat.A32, mat.A33 = rib_dir.z, span_dir.z, third_dir.z
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
        except:
            pass
    return unified

def create_cut_result(body, params, doc=None):
    """Returns (cut_shape, rib_center_lines, list_of_valid_rib_solids)"""
    if doc is None:
        doc = FreeCAD.ActiveDocument
    raw_shape = body.Shape if hasattr(body, "Shape") else body
    bb = raw_shape.BoundBox

    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_center_lines = lines1 + lines2

    rib_faces1 = create_rib_faces(lines1, params.plane_normal, params.rib_width)
    rib_faces2 = create_rib_faces(lines2, params.plane_normal, params.rib_width)

    rib_solids1 = extrude_rib_faces_to_solids(rib_faces1, params.plane_normal, bb)
    rib_solids2 = extrude_rib_faces_to_solids(rib_faces2, params.plane_normal, bb)
    all_rib_solids = rib_solids1 + rib_solids2

    # ------------------------------------------------------------
    # FILTER: keep only ribs that significantly intersect the wing
    # ------------------------------------------------------------
    valid_ribs = []
    for rib in all_rib_solids:
        try:
            # Intersect the rib solid with the wing
            common = rib.common(raw_shape)
            # If the common shape has volume less than 1 mm³, discard it
            if common is None or common.Volume < 1.0:
                continue
            valid_ribs.append(rib)
        except Exception:
            continue

    if not valid_ribs:
        raise ValueError("No rib solids significantly intersect the wing.")

    # Fuse the valid ribs for cutting
    def fuse_list(lst):
        if not lst:
            return None
        f = lst[0]
        for s in lst[1:]:
            try:
                f = f.fuse(s)
            except Exception:
                pass
        return f

    fused_ribs = fuse_list(valid_ribs)
    cut_result = raw_shape.cut(fused_ribs)

    return cut_result, all_center_lines, valid_ribs

def generate_final_solid(body, params, doc=None, add_bridges=True):
    """Returns the final solid (cut + bridges)"""
    cut_result, center_lines, rib_solids = create_cut_result(body, params, doc)
    if add_bridges:
        result = add_bridges_along_ribs(cut_result, center_lines,
                                        params.rib_width, params.plane_normal,
                                        bridge_size=params.rib_width * 10)
    else:
        result = cut_result
    if doc:
        obj = doc.addObject("Part::Feature", "WingWithInfill")
        obj.Shape = result
        doc.recompute()
    return result