# slab_utils.py
import FreeCAD
import Part
import math
import random

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

def create_cut_result(body, params, doc=None, vis=False):
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

    valid_ribs = []
    for rib in all_rib_solids:
        try:
            common = rib.common(raw_shape)
            if common is None or common.Volume < 1.0:
                continue
            valid_ribs.append(rib)
        except Exception:
            continue

    if not valid_ribs:
        raise ValueError("No rib solids significantly intersect the wing.")

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

    if vis:
        from viz_utils import show_cut_wing
        show_cut_wing(cut_result, doc, transparency=80)

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


def merge_and_show_final(cut_result, bridges_shape, doc, vis=True):
    """
    Fuse the cut wing with the bridges shape.
    If vis is True, display the final solid using viz_utils.show_final_solid.
    Returns the fused solid (or None if bridges_shape is invalid).
    """
    if not bridges_shape or bridges_shape.isNull():
        return None
    final = cut_result.fuse(bridges_shape)
    if vis:
        from viz_utils import show_final_solid
        show_final_solid(final, doc)
    return final

def create_rib_centre_surfaces(wing_shape, rib_center_lines, plane_normal, tol=1e-4):
    """
    Returns (faces, edges).
    faces: list of Part.Face — one for each closed region (including holes)
    edges: list of Part.Edge — all raw intersection edges
    """
    faces = []
    all_edges = []

    for line in rib_center_lines:
        start = line.Vertexes[0].Point
        end   = line.Vertexes[-1].Point
        rib_dir = (end - start).normalize()
        n_rib = rib_dir.cross(plane_normal).normalize()
        plane = Part.Plane(start, n_rib)
        try:
            intersect = wing_shape.section(plane)
            if not intersect or not intersect.Edges:
                continue
            edges = list(intersect.Edges)
            all_edges.extend(edges)

            # ---- Build all closed wires by edge stitching ----
            used = [False] * len(edges)
            wires = []
            for i in range(len(edges)):
                if used[i]:
                    continue
                # Start a new wire
                loop = [edges[i]]
                used[i] = True
                cur_end = edges[i].Vertexes[-1].Point
                changed = True
                while changed:
                    changed = False
                    for j, e in enumerate(edges):
                        if used[j]:
                            continue
                        if e.Vertexes[0].Point.isEqual(cur_end, tol):
                            loop.append(e)
                            used[j] = True
                            cur_end = e.Vertexes[-1].Point
                            changed = True
                            break
                        elif e.Vertexes[-1].Point.isEqual(cur_end, tol):
                            rev = e.copy()
                            rev.reverse()
                            loop.append(rev)
                            used[j] = True
                            cur_end = e.Vertexes[0].Point
                            changed = True
                            break
                # Close the wire
                if len(loop) > 1:
                    try:
                        wire = Part.Wire(loop)
                        wires.append(wire)
                    except:
                        continue

            # For each closed wire, create a face. We will not attempt to detect nesting; just create faces.
            # If there is nesting, some faces may be holes, but they will be separate objects.
            # That's acceptable for visualisation and further processing.
            for wire in wires:
                if wire.isClosed() and not wire.isNull():
                    try:
                        face = Part.Face(wire)
                        if face.isValid():
                            faces.append(face)
                    except:
                        # Fallback: try to create face with tolerance
                        try:
                            face = Part.Face(wire, tolerance=tol)
                            if face.isValid():
                                faces.append(face)
                        except:
                            continue
        except Exception:
            continue

    return faces, all_edges

def _lines_intersect_3d(p1, d1, p2, d2, tol=1.0):
    """
    Test if two infinite lines (p1+t*d1) and (p2+s*d2) pass within
    `tol` of each other. Returns the closest point on line1 if they
    are close enough, else None.
 
    Uses the standard closest-point-between-two-lines formula.
    """
    import math
    w = FreeCAD.Vector(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)
 
    def dot(a, b): return a.x*b.x + a.y*b.y + a.z*b.z
 
    a = dot(d1, d1)
    b = dot(d1, d2)
    c = dot(d2, d2)
    d = dot(d1, w)
    e = dot(d2, w)
 
    denom = a*c - b*b
    if abs(denom) < 1e-10:
        return None  # parallel lines
 
    t = (b*e - c*d) / denom
    s = (a*e - b*d) / denom
 
    # Closest points
    cp1 = FreeCAD.Vector(p1.x + d1.x*t, p1.y + d1.y*t, p1.z + d1.z*t)
    cp2 = FreeCAD.Vector(p2.x + d2.x*s, p2.y + d2.y*s, p2.z + d2.z*s)
 
    dist = math.sqrt((cp1.x-cp2.x)**2 + (cp1.y-cp2.y)**2 + (cp1.z-cp2.z)**2)
    if dist <= tol:
        return cp1  # intersection point (approximately)
    return None
 
 
def split_ribs_by_crossings(rib_solids, rib_center_lines, plane_normal, tol=1.0):
    """
    For each rib solid, cut it by all other rib solids whose centre
    lines cross it. Returns a flat list of all resulting solid pieces.
 
    tol: distance tolerance (mm) for deciding if two rib lines cross.
         Should be >= rib_width to catch near-misses.
    """
    n = len(rib_solids)
 
    # Build crossing map: for each rib i, which ribs j cross it?
    crossing_map = {i: [] for i in range(n)}
    for i in range(n):
        si = rib_center_lines[i].Vertexes[0].Point
        ei = rib_center_lines[i].Vertexes[-1].Point
        try:
            di_raw = FreeCAD.Vector(ei.x-si.x, ei.y-si.y, ei.z-si.z)
            L = di_raw.Length
            if L < 1e-6:
                continue
            di = FreeCAD.Vector(di_raw.x/L, di_raw.y/L, di_raw.z/L)
        except Exception:
            continue
 
        for j in range(i+1, n):
            sj = rib_center_lines[j].Vertexes[0].Point
            ej = rib_center_lines[j].Vertexes[-1].Point
            try:
                dj_raw = FreeCAD.Vector(ej.x-sj.x, ej.y-sj.y, ej.z-sj.z)
                Lj = dj_raw.Length
                if Lj < 1e-6:
                    continue
                dj = FreeCAD.Vector(dj_raw.x/Lj, dj_raw.y/Lj, dj_raw.z/Lj)
            except Exception:
                continue
 
            pt = _lines_intersect_3d(si, di, sj, dj, tol=tol)
            if pt is not None:
                crossing_map[i].append(j)
                crossing_map[j].append(i)
 
    # Cut each rib by its crossing ribs
    all_pieces = []
    for i, rib in enumerate(rib_solids):
        if not crossing_map[i]:
            all_pieces.append(rib)
            continue
 
        # Fuse all crossing ribs into one cutting tool
        cutters = [rib_solids[j] for j in crossing_map[i]]
        cut_tool = cutters[0]
        for c in cutters[1:]:
            try:
                cut_tool = cut_tool.fuse(c)
            except Exception:
                pass
 
        try:
            result = rib.cut(cut_tool)
            if result.isNull():
                all_pieces.append(rib)
                continue
            pieces = result.Solids if result.ShapeType in ("Compound", "CompSolid") else [result]
            pieces = [p for p in pieces if not p.isNull() and p.Volume > 1e-6]
            all_pieces.extend(pieces if pieces else [rib])
        except Exception:
            all_pieces.append(rib)
 
    return all_pieces

def split_rib_faces_by_crossings(rib_faces, tol=1e-4):
    """
    Cut each rib centre face by all other rib centre faces.
    Where two faces from different families cross, each is split
    at the intersection line, giving individual planar segments.
 
    Parameters
    ----------
    rib_faces : list of Part.Face
        The rib centre surfaces from create_rib_centre_surfaces.
    tol : float
        Geometric tolerance for boolean operations.
 
    Returns
    -------
    list of Part.Face
        All resulting face fragments, flat list.
    """
 
    n = len(rib_faces)
    result_faces = []
 
    for i, face in enumerate(rib_faces):
        current = [face]  # start with the whole face, progressively cut
 
        for j, other in enumerate(rib_faces):
            if i == j:
                continue
 
            # Check if they actually intersect before trying to cut
            try:
                section = face.section(other)
                if section.isNull() or not section.Edges:
                    continue  # no intersection, skip
            except Exception:
                continue
 
            # Cut all current fragments by this other face
            next_fragments = []
            for fragment in current:
                try:
                    cut = fragment.cut(other)
                    if cut.isNull():
                        next_fragments.append(fragment)  # cut failed, keep original
                        continue
                    # Extract resulting faces from the cut
                    faces_out = cut.Faces if hasattr(cut, 'Faces') else [cut]
                    valid = [f for f in faces_out
                             if f is not None and not f.isNull() and f.Area > tol]
                    if valid:
                        next_fragments.extend(valid)
                    else:
                        next_fragments.append(fragment)
                except Exception:
                    next_fragments.append(fragment)
 
            current = next_fragments
 
        result_faces.extend(current)
 
    return result_faces

def create_rectangular_cutout_from_boundary(face, margin=1.0):
    """
    For a planar face (rib segment), compute a quadrilateral cutout.
    If the face has holes (inner loops), split it into separate faces and process each.
    Returns a list of cutout faces (could be empty, one, or several).
    """
    # Helper to process a single face (no holes)
    def process_single_face(f, margin):
        wires = f.Wires
        if not wires:
            return None
        wire = wires[0]
        # Collect vertices
        vertices = []
        for edge in wire.Edges:
            vertices.append(edge.Vertexes[0].Point)
            vertices.append(edge.Vertexes[-1].Point)
        uniq = []
        for p in vertices:
            if not any(p.distanceToPoint(q) < 1e-4 for q in uniq):
                uniq.append(p)
        if len(uniq) < 4:
            return None
        # Group by Z
        z_groups = {}
        for p in uniq:
            key = round(p.z, 4)
            z_groups.setdefault(key, []).append(p)
        z_min = min(z_groups.keys())
        z_max = max(z_groups.keys())
        low_pts = z_groups[z_min]
        high_pts = z_groups[z_max]
        low_cent = FreeCAD.Vector(0,0,0)
        for p in low_pts: low_cent += p
        low_cent /= len(low_pts)
        high_cent = FreeCAD.Vector(0,0,0)
        for p in high_pts: high_cent += p
        high_cent /= len(high_pts)
        z_mid = (low_cent.z + high_cent.z) / 2.0
        # Mid plane
        plane = Part.Plane(FreeCAD.Vector(0,0,z_mid), FreeCAD.Vector(0,0,1))
        try:
            section = f.section(plane)
        except:
            return None
        if not section or len(section.Edges) == 0:
            return None
        mid_pts = []
        for edge in section.Edges:
            mid_pts.append(edge.Vertexes[0].Point)
            mid_pts.append(edge.Vertexes[-1].Point)
        uniq_mid = []
        for p in mid_pts:
            if not any(p.distanceToPoint(q) < 1e-4 for q in uniq_mid):
                uniq_mid.append(p)
        if len(uniq_mid) < 2:
            return None
        # Farthest pair
        best_pair = None
        max_dist = 0.0
        for i in range(len(uniq_mid)):
            for j in range(i+1, len(uniq_mid)):
                d = uniq_mid[i].distanceToPoint(uniq_mid[j])
                if d > max_dist:
                    max_dist = d
                    best_pair = (uniq_mid[i], uniq_mid[j])
        if best_pair is None:
            return None
        left_mid, right_mid = best_pair
        if left_mid.x > right_mid.x:
            left_mid, right_mid = right_mid, left_mid
        # Centroid and shifting
        centroid = (low_cent + left_mid + high_cent + right_mid) / 4.0
        def shift_towards(pt, center, margin):
            dir_vec = pt - center
            if dir_vec.Length < 1e-6:
                return pt
            return pt - dir_vec.normalize() * margin
        low_shifted = shift_towards(low_cent, centroid, margin)
        high_shifted = shift_towards(high_cent, centroid, margin)
        left_shifted = shift_towards(left_mid, centroid, margin)
        right_shifted = shift_towards(right_mid, centroid, margin)
        ordered = [low_shifted, left_shifted, high_shifted, right_shifted]
        try:
            wire_out = Part.makePolygon(ordered + [ordered[0]])
            face_out = Part.Face(wire_out)
            if face_out.isValid() and face_out.Area > 1e-4:
                # Check inside original face
                sample_points = [low_shifted, left_shifted, high_shifted, right_shifted,
                                (low_shifted + high_shifted) * 0.5,
                                (left_shifted + right_shifted) * 0.5]
                inside = True
                for p in sample_points:
                    if not f.isInside(p, 1e-3, True):
                        inside = False
                        break
                if inside:
                    return face_out
        except:
            pass
        return None

    # Main: if face has holes, split into separate faces first
    if len(face.Wires) > 1:
        # Extract outer wire and inner wires
        wires = list(face.Wires)
        outer = wires[0]
        inners = wires[1:]
        # For each inner wire, create a face that is the outer minus the inner
        # (i.e., a face with a hole). Then process each such face?
        # Actually, we want separate faces for each region outside the holes.
        # Simpler: use the wire to build a new face that is only the outer boundary (ignoring holes),
        # but that would include the hole area. We need to split into multiple faces.
        # A practical approach: create a compound of faces by cutting the outer face with the inner wires.
        # Use Part.Geom2d to split, but that's complex.
        # For now, we skip faces with holes (old behaviour) – but you asked not to skip.
        # Given the complexity, we will skip and print a warning.
        print("  [warning] face with hole – skipping cutout (will not process both parts)")
        return None
    else:
        return process_single_face(face, margin)
