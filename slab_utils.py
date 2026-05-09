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
    faces: list of Part.Face — one per rib where a surface could be built
    edges: list of Part.Edge — all raw intersection edges (always returned)

    Strategy per rib:
      1. Intersect rib plane with wing solid via section()
      2. Collect all edges and add to all_edges
      3. Try to stitch edges into the longest connected wire
      4. Try Part.Face (requires closed planar wire)
      5. Fallback: Part.makeFilledFace (handles open/non-planar)
    """
    faces = []
    all_edges = []

    for line in rib_center_lines:
        start   = line.Vertexes[0].Point
        end     = line.Vertexes[-1].Point
        rib_dir = (end - start).normalize()
        n_rib   = rib_dir.cross(plane_normal).normalize()
        plane   = Part.Plane(start, n_rib)

        try:
            intersect = wing_shape.section(plane)
            if not intersect or not intersect.Edges:
                continue

            edges = list(intersect.Edges)
            all_edges.extend(edges)

            # ---- Stitch edges into best connected loop ----
            used = [False] * len(edges)
            best_loop = []

            for i in range(len(edges)):
                if used[i]:
                    continue
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
                if len(loop) > len(best_loop):
                    best_loop = loop

            if not best_loop:
                continue

            # ---- Build wire ----
            try:
                wire = Part.Wire(best_loop)
            except Exception:
                continue

            # ---- Try Part.Face (closed planar wire) ----
            face = None
            if wire.isClosed():
                try:
                    f = Part.Face(wire)
                    if f.isValid():
                        face = f
                except Exception:
                    pass

            # ---- Fallback: Part.makeFilledFace ----
            if face is None:
                try:
                    f = Part.makeFilledFace(wire.Edges)
                    if f.isValid():
                        face = f
                except Exception:
                    pass

            if face is not None:
                faces.append(face)

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

def add_ellipse_holes_to_faces(faces, margin=0.5):
    """
    For each face, compute the largest inscribed ellipse with a given
    margin (minimum distance from the edge) and return a list of (face, ellipse_face)
    pairs, or directly cut the face and return the holed face.
    """
    from math import cos, sin, pi
    
    def distance_to_boundary(p, boundary_points):
        """Minimum distance from point p (2D Vector) to any boundary edge."""
        best = float('inf')
        for i in range(len(boundary_points)):
            a = boundary_points[i]
            b = boundary_points[(i+1)%len(boundary_points)]
            # Distance from point to segment ab
            ab = b - a
            t = (p - a).dot(ab) / ab.dot(ab)
            if t < 0:
                d = (p - a).Length
            elif t > 1:
                d = (p - b).Length
            else:
                proj = a + ab * t
                d = (p - proj).Length
            if d < best:
                best = d
        return best
    
    holed_faces = []
    for face in faces:
        # Get outer wire (assume face has one wire)
        wires = face.Wires
        if not wires:
            continue
        wire = wires[0]
        # Extract boundary vertices in 2D (local coordinates)
        # We'll work in the plane of the face: using face.Surface.parameter(p)
        # But simpler: project points to a 2D coordinate system on the face.
        # Use the face's own parameterisation.
        # For a planar face, we can use face.Surface.param(p) which returns (u,v).
        # However, the distance in 3D is proportional to distance in UV if the face is isometric.
        # For simplicity, we'll sample the 3D points and project to a local 2D coordinate system.
        # Get a normal and two orthogonal axes.
        centre = face.CenterOfMass
        normal = face.normalAt(0.5,0.5).normalize()
        # Choose X and Y axes arbitrarily in the plane
        # Pick a vector not parallel to normal
        if abs(normal.x) < 0.9:
            x_axis = FreeCAD.Vector(1,0,0).cross(normal).normalize()
        else:
            x_axis = FreeCAD.Vector(0,1,0).cross(normal).normalize()
        y_axis = normal.cross(x_axis).normalize()
        # Convert 3D points to 2D
        boundary_2d = []
        for edge in wire.Edges:
            pts = edge.discretize(20)  # sample each edge
            for p in pts:
                dx = (p - centre).dot(x_axis)
                dy = (p - centre).dot(y_axis)
                boundary_2d.append(FreeCAD.Vector(dx, dy, 0))
        # Remove duplicates
        uniq = []
        for p in boundary_2d:
            if not any(p.distanceToPoint(q) < 1e-6 for q in uniq):
                uniq.append(p)
        boundary_2d = uniq
        if len(boundary_2d) < 3:
            continue
        
        # Sample candidate points inside (grid)
        bbox_min = FreeCAD.Vector(min(p.x for p in boundary_2d), min(p.y for p in boundary_2d), 0)
        bbox_max = FreeCAD.Vector(max(p.x for p in boundary_2d), max(p.y for p in boundary_2d), 0)
        step = 2.0   # mm – adjust for resolution
        best_center = None
        best_radius = 0.0
        x = bbox_min.x
        while x <= bbox_max.x:
            y = bbox_min.y
            while y <= bbox_max.y:
                p = FreeCAD.Vector(x, y, 0)
                # Quick bounding box test
                if p.x < bbox_min.x or p.x > bbox_max.x or p.y < bbox_min.y or p.y > bbox_max.y:
                    y += step
                    continue
                dist = distance_to_boundary(p, boundary_2d)
                if dist > best_radius:
                    best_radius = dist
                    best_center = p
                y += step
            x += step
        
        if best_center is None or best_radius < margin:
            continue  # no space for hole
        
        # Reduce radius by margin
        radius = best_radius - margin
        if radius < 0.5:
            continue
        
        # Create an ellipse (circle for simplicity) in the plane
        # For ellipse, we can use the sampled distances along two perpendicular directions.
        # Here we'll just make a circle (special case of ellipse).
        # Convert back to 3D
        center_3d = centre + x_axis * best_center.x + y_axis * best_center.y
        # Create a face for the ellipse
        circle = Part.makeCircle(radius, center_3d, normal)
        ellipse_wire = Part.Wire(circle)
        ellipse_face = Part.Face(ellipse_wire)
        # Cut the original face with the ellipse
        try:
            holed = face.cut(ellipse_face)
            if not holed.isNull():
                holed_faces.append(holed)
            else:
                holed_faces.append(face)
        except:
            holed_faces.append(face)
    return holed_faces

def is_point_inside_polygon(p, poly):
    """Ray casting algorithm for 2D point."""
    x, y = p.x, p.y
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i].x, poly[i].y
        x2, y2 = poly[(i+1)%n].x, poly[(i+1)%n].y
        if ((y1 > y) != (y2 > y)) and (x < (x2-x1)*(y-y1)/(y2-y1)+x1):
            inside = not inside
    return inside

def ellipse_distance_to_boundary(cx, cy, a, b, ang, boundary_points, margin, n_samples=36):
    """
    Compute min distance from ellipse boundary to the polygon edges,
    plus a penalty if ellipse centre is outside.
    """
    # First, ensure centre is inside polygon (quick guard)
    if not is_point_inside_polygon(FreeCAD.Vector(cx,cy,0), boundary_points):
        return -1.0  # invalid
    
    # Sample points on ellipse
    best = float('inf')
    for i in range(n_samples):
        t = 2*math.pi * i / n_samples
        x = cx + a * math.cos(t) * math.cos(ang) - b * math.sin(t) * math.sin(ang)
        y = cy + a * math.cos(t) * math.sin(ang) + b * math.sin(t) * math.cos(ang)
        pt = FreeCAD.Vector(x, y, 0)
        # Distance to nearest boundary edge
        dist = min(point_edge_distance(pt, boundary_points) for _ in range(1))
        # Actually compute properly:
        min_d = float('inf')
        n = len(boundary_points)
        for j in range(n):
            p1 = boundary_points[j]
            p2 = boundary_points[(j+1)%n]
            # Distance from pt to segment p1-p2
            ab = p2 - p1
            t_param = (pt - p1).dot(ab) / ab.dot(ab)
            if t_param < 0:
                d = (pt - p1).Length
            elif t_param > 1:
                d = (pt - p2).Length
            else:
                proj = p1 + ab * t_param
                d = (pt - proj).Length
            if d < min_d:
                min_d = d
        if min_d < best:
            best = min_d
        if best < 0:
            break
    # We want distance >= margin, so we minimise negative margin distance
    # objective = -min_distance (so larger distance better)
    # but we need to ensure margin satisfied.
    return best  # min distance from ellipse to boundary

def point_edge_distance(pt, poly):
    """Minimum distance from pt to any edge of polygon."""
    min_d = float('inf')
    n = len(poly)
    for i in range(n):
        p1 = poly[i]
        p2 = poly[(i+1)%n]
        ab = p2 - p1
        t = (pt - p1).dot(ab) / ab.dot(ab)
        if t < 0:
            d = (pt - p1).Length
        elif t > 1:
            d = (pt - p2).Length
        else:
            proj = p1 + ab * t
            d = (pt - proj).Length
        if d < min_d:
            min_d = d
    return min_d

def find_largest_ellipse(boundary_points, margin=1.0, steps=20, restarts=5):
    """
    Returns (cx, cy, a, b, ang) of largest ellipse with distance >= margin.
    """
    # Bounding box of polygon
    xs = [p.x for p in boundary_points]
    ys = [p.y for p in boundary_points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    
    best_obj = -1e9
    best_params = None
    
    for _ in range(restarts):
        # Random initial centre inside polygon
        while True:
            cx = random.uniform(xmin, xmax)
            cy = random.uniform(ymin, ymax)
            if is_point_inside_polygon(FreeCAD.Vector(cx,cy,0), boundary_points):
                break
        # Initial axes as fraction of bounding box
        a_init = (xmax - xmin) * random.uniform(0.1, 0.4)
        b_init = (ymax - ymin) * random.uniform(0.1, 0.4)
        ang_init = random.uniform(0, math.pi)
        
        # Local optimisation: simple hill‑climbing
        # Parameters: [cx, cy, a, b, ang]
        params = [cx, cy, a_init, b_init, ang_init]
        for iteration in range(steps):
            # Try small perturbations in each direction
            best_local = None
            best_val = ellipse_distance_to_boundary(params[0], params[1], params[2], params[3], params[4], boundary_points, margin)
            for dp in [0.5, -0.5, 1.0, -1.0, 2.0, -2.0]:  # varying step sizes
                for i in range(5):
                    new_params = params[:]
                    new_params[i] += dp
                    # Clamp a,b to positive
                    if new_params[2] < 0.5: new_params[2] = 0.5
                    if new_params[3] < 0.5: new_params[3] = 0.5
                    # Clamp centre to stay inside bounding box (approx)
                    new_params[0] = max(xmin+0.5, min(xmax-0.5, new_params[0]))
                    new_params[1] = max(ymin+0.5, min(ymax-0.5, new_params[1]))
                    val = ellipse_distance_to_boundary(new_params[0], new_params[1], new_params[2], new_params[3], new_params[4], boundary_points, margin)
                    if val > best_val:
                        best_val = val
                        best_local = new_params[:]
            if best_local:
                params = best_local
            else:
                break
        # After local search, compute area
        if best_val > best_obj and best_val >= margin:
            best_obj = best_val
            best_params = params[:]
    if best_params is None:
        return None
    # Return also the area = pi*a*b
    area = math.pi * best_params[2] * best_params[3]
    return best_params + [area]
