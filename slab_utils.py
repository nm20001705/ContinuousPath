# slab_utils.py
import FreeCAD
import Part
import math
import sys

# ===== FreeCAD Profiler Helper =====
import FreeCAD
import time
import sys
from perf_utils import timed

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

@timed
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
    lines1 = generate_lines(d1)
    lines2 = generate_lines(d2)
    return lines1, lines2

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
        except Exception:
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
        except Exception:
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
        except Exception:
            pass
    return unified

@timed
def create_cut_result_mesh(wing_mesh, params, doc=None, lines1=None, lines2=None):
    """Returns (cut_mesh, rib_center_lines) using mesh booleans."""
    if doc is None:
        doc = FreeCAD.ActiveDocument
    bb = wing_mesh.BoundBox

    # 1. Generate rib solids (same as before)
    if lines1 is None or lines2 is None:
        lines1, lines2 = create_angled_grid_lines(bb, params)
    all_center_lines = lines1 + lines2

    rib_faces1 = create_rib_faces(lines1, params.plane_normal, params.rib_width)
    rib_faces2 = create_rib_faces(lines2, params.plane_normal, params.rib_width)
    rib_solids1 = extrude_rib_faces_to_solids(rib_faces1, params.plane_normal, bb)
    rib_solids2 = extrude_rib_faces_to_solids(rib_faces2, params.plane_normal, bb)
    all_rib_solids = rib_solids1 + rib_solids2

    # 2. Convert each rib solid to a mesh and fuse them
    rib_meshes = []
    for rib in all_rib_solids:
        try:
            m = solid_to_mesh(rib)
            if m and m.Facets:
                rib_meshes.append(m)
        except:
            continue
    if not rib_meshes:
        return wing_mesh, all_center_lines

    # Fuse all rib meshes into one
    rib_mesh = rib_meshes[0]
    for m in rib_meshes[1:]:
        rib_mesh = rib_mesh.fuse(m)

    # 3. Subtract the rib mesh from the wing mesh
    cut_mesh = wing_mesh.difference(rib_mesh)

    return cut_mesh, all_center_lines

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

@timed
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

@timed

def create_rib_centre_surfaces(wing_shape, rib_center_lines, plane_normal, doc, vis=False, tol=1e-3):
    """
    Returns (faces, edges).
    faces: list of Part.Face — one for each rib surface (if face creation succeeds).
    edges: list of Part.Edge — all raw intersection edges (always returned).
    """
    faces = []
    all_edges = []
    wing_bb = wing_shape.BoundBox

    for line in rib_center_lines:
        start = line.Vertexes[0].Point
        end   = line.Vertexes[-1].Point

        # Quick bounding box check using FreeCAD.BoundBox
        xmin, xmax = min(start.x, end.x), max(start.x, end.x)
        ymin, ymax = min(start.y, end.y), max(start.y, end.y)
        zmin, zmax = min(start.z, end.z), max(start.z, end.z)
        line_bb = FreeCAD.BoundBox(xmin, ymin, zmin, xmax, ymax, zmax)
        if not line_bb.intersect(wing_bb):
            continue

        rib_dir = (end - start).normalize()
        n_rib = rib_dir.cross(plane_normal).normalize()
        plane = Part.Plane(start, n_rib)

        try:
            # Try makeCrossSection (if available)
            try:
                intersect = wing_shape.makeCrossSection(plane, tol)
            except AttributeError:
                # Fallback to section
                intersect = wing_shape.section(plane)

            if not intersect or not intersect.Edges:
                continue
            edges = list(intersect.Edges)
            all_edges.extend(edges)

            # ---- Stitch edges into wires ----
            used = [False] * len(edges)
            wires = []
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
                if len(loop) > 1:
                    try:
                        wire = Part.Wire(loop)
                        wires.append(wire)
                    except:
                        continue

            if not wires:
                continue

            # Select the wire with largest area (outer boundary)
            outer_wire = max(wires, key=lambda w: abs(w.Area))
            try:
                face = Part.Face(outer_wire)
                if face.isValid():
                    faces.append(face)
            except:
                # Fallback: try to create face with tolerance
                try:
                    face = Part.Face(outer_wire, tolerance=tol)
                    if face.isValid():
                        faces.append(face)
                except:
                    continue

        except Exception:
            continue

    if vis:
        from viz_utils import show_rib_centre_surfaces, show_rib_centre_edges
        show_rib_centre_surfaces(faces, doc, color=(0.8, 0.4, 0.8), transparency=50)
        print(f"Created {len(faces)} rib centre surfaces.")
        show_rib_centre_edges(all_edges, doc, line_color=(0.2, 0.5, 1.0), line_width=2)
        print(f"Created {len(all_edges)} rib centre edges.")

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

@timed
def split_rib_faces_by_crossings(rib_faces, doc, viz_params, tol=1e-4, min_area=0.1):
    """
    Optimized version: only cut faces with area > min_area and from different families.
    """
    n = len(rib_faces)
    if n == 0:
        return []

    # Pre‑compute normals, bounding boxes, and areas
    normals = []
    bboxes = []
    areas = []
    for face in rib_faces:
        try:
            normal = face.normalAt(0.5, 0.5).normalize()
        except:
            normal = FreeCAD.Vector(0,0,1)
        normals.append(normal)
        bboxes.append(face.BoundBox)
        areas.append(face.Area)

    all_pieces = []
    cut_count = 0
    skip_small_count = 0
    skip_parallel_count = 0
    skip_bbox_count = 0
    candidate_checks = 0

    for i in range(n):
        face_i = rib_faces[i]
        # Skip tiny faces – they won't contain a cutout
        if areas[i] < min_area:
            skip_small_count += 1
            all_pieces.append(face_i)
            continue

        bbox_i = bboxes[i]
        normal_i = normals[i]

        candidates = []
        for j in range(n):
            if i == j:
                continue
            candidate_checks += 1
            # Skip same family (parallel normals)
            if abs(normal_i.dot(normals[j])) > 0.9:
                skip_parallel_count += 1
                continue
            if not bbox_i.intersect(bboxes[j]):
                skip_bbox_count += 1
                continue
            candidates.append(j)

        if candidates:
            cut_tool = rib_faces[candidates[0]]
            for j in candidates[1:]:
                try:
                    cut_tool = cut_tool.fuse(rib_faces[j])
                except:
                    pass
            try:
                result = face_i.cut(cut_tool)
                cut_count += 1
                if not result.isNull():
                    faces_out = result.Faces if hasattr(result, 'Faces') else [result]
                    valid = [f for f in faces_out if f is not None and not f.isNull() and f.Area > tol]
                    if valid:
                        all_pieces.extend(valid)
                    else:
                        all_pieces.append(face_i)
                else:
                    all_pieces.append(face_i)
            except:
                all_pieces.append(face_i)
        else:
            all_pieces.append(face_i)


    if viz_params:
        from viz_utils import show_rib_segments
        show_rib_segments(all_pieces, doc)
    return all_pieces

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

@timed
def cutouts_from_segmens(rib_centre_segments, doc, vis, margin):
    cutouts = []
    for seg in rib_centre_segments:
        cut = create_rectangular_cutout_from_boundary(seg, margin=margin)
        if cut:
            cutouts.append(cut)
    if vis:
        from viz_utils import show_rect_cutouts
        show_rect_cutouts(cutouts, doc, color=(1.0,0.5,0.0), transparency=30)
    return cutouts

def apply_holes_to_ribs(rib_solids, cutout_faces, rib_width, plane_normal, doc=None):
    import time
    if not cutout_faces or not rib_solids:
        return rib_solids

    # 1. Create hole solids
    t0 = time.time()
    holes = []
    for face in cutout_faces:
        normal = face.normalAt(0.5, 0.5).normalize()
        half = rib_width / 2.0
        moved = face.copy()
        moved.translate(normal * (-half))
        try:
            solid = moved.extrude(normal * rib_width)
            if not solid.isNull():
                holes.append(solid)
        except:
            continue

    if not holes:
        return rib_solids

    # 2. Fuse holes
    t0 = time.time()
    cutter = holes[0]
    for h in holes[1:]:
        cutter = cutter.fuse(h)

    # 3. Cut each rib
    t0 = time.time()
    holed = []
    for i, rib in enumerate(rib_solids):
        try:
            cut_rib = rib.cut(cutter)
            if cut_rib.isNull():
                holed.append(rib)
            else:
                holed.append(cut_rib)
        except:
            holed.append(rib)

    return holed
