# point_utils.py
import FreeCAD
import Part
import math

def _get_wing_wires_at_z(wing_shape, z):
    """
    Return a list of closed wires representing the cross-section of the wing
    at height z using the slice() method.
    """
    try:
        # slice(direction, distance) returns a compound or list of compounds
        sliced = wing_shape.slice(FreeCAD.Vector(0, 0, 1), z)
        if not sliced:
            return []
        # If result is a list, iterate; otherwise treat as a compound
        if isinstance(sliced, list):
            compounds = sliced
        else:
            compounds = [sliced]
        wires = []
        for comp in compounds:
            if hasattr(comp, 'Edges') and len(comp.Edges) > 0:
                # Try to build wires from edges
                try:
                    # First, try to create a single wire from all edges (if connected)
                    wire = Part.Wire(comp.Edges)
                    wires.append(wire)
                except:
                    # Edges may belong to multiple closed loops – we need to sort them
                    # For simplicity, we create a compound of all edges – later we'll treat them as a single wire
                    # but better: sew edges into wires
                    from collections import defaultdict
                    # create a list of edges
                    edges = list(comp.Edges)
                    used = [False]*len(edges)
                    for i in range(len(edges)):
                        if used[i]:
                            continue
                        wire_edges = [edges[i]]
                        used[i] = True
                        cur_end = edges[i].Vertexes[-1].Point
                        changed = True
                        while changed:
                            changed = False
                            for j, e in enumerate(edges):
                                if used[j]:
                                    continue
                                if e.Vertexes[0].Point.isEqual(cur_end, 1e-6):
                                    wire_edges.append(e)
                                    used[j] = True
                                    cur_end = e.Vertexes[-1].Point
                                    changed = True
                                    break
                                elif e.Vertexes[-1].Point.isEqual(cur_end, 1e-6):
                                    # reverse edge
                                    rev = Part.Edge(e.Curve, e.Vertexes[-1].Point, e.Vertexes[0].Point)
                                    wire_edges.append(rev)
                                    used[j] = True
                                    cur_end = e.Vertexes[0].Point
                                    changed = True
                                    break
                        if len(wire_edges) > 2:
                            try:
                                w = Part.Wire(wire_edges)
                                wires.append(w)
                            except:
                                pass
        return wires
    except Exception as e:
        print(f"Error getting wires at z={z}: {e}")
        return []

def _rib_plane_intersect_z_plane(rib_center_line, plane_normal, z):
    """
    Compute intersection of rib slab plane with z=const.
    Returns (point_on_line, line_dir) or None.
    """
    start = rib_center_line.Vertexes[0].Point
    end   = rib_center_line.Vertexes[-1].Point
    d = (end - start).normalize()
    # Rib plane normal = d × plane_normal
    n_rib = d.cross(plane_normal).normalize()
    Z = FreeCAD.Vector(0,0,1)
    line_dir = n_rib.cross(Z).normalize()
    if line_dir.Length < 1e-6:
        return None
    c = n_rib.dot(start)
    rhs = c - n_rib.z * z
    if abs(n_rib.x) > abs(n_rib.y):
        if abs(n_rib.x) < 1e-6:
            return None
        x = rhs / n_rib.x
        y = 0.0
    else:
        if abs(n_rib.y) < 1e-6:
            return None
        x = 0.0
        y = rhs / n_rib.y
    point_on_line = FreeCAD.Vector(x, y, z)
    return point_on_line, line_dir

def _clip_line_to_wing_wires(point_on_line, line_dir, wires, max_len=2000):
    if not wires:
        return None
    p = point_on_line
    d = line_dir
    seg_start = p - d * max_len
    seg_end   = p + d * max_len
    test_edge = Part.makeLine(seg_start, seg_end)
    t_vals = []
    for wire in wires:
        for edge in wire.Edges:
            try:
                dist, pts, _ = test_edge.distToShape(edge)
                if dist < 0.1:
                    for pt_pair in pts:
                        t = (pt_pair[0] - p).dot(d)
                        t_vals.append(t)
            except:
                continue
    if len(t_vals) < 2:
        return None
    t_vals.sort()
    best_seg = None
    best_len = 0.0
    for i in range(len(t_vals)-1):
        t1, t2 = t_vals[i], t_vals[i+1]
        t_mid = (t1+t2)/2.0
        mid_pt = p + d * t_mid
        inside = False
        for wire in wires:
            try:
                face = Part.Face(wire)
                if face.isInside(mid_pt, 1e-6, True):
                    inside = True
                    break
            except:
                if wire.BoundBox.isInside(mid_pt):
                    inside = True
                    break
        if inside:
            length = t2 - t1
            if length > best_len:
                best_len = length
                best_seg = (p + d * t1, p + d * t2)
    return best_seg

def collect_rib_midpoints(wing_shape, rib_center_lines, plane_normal, z_min, z_max, z_step):
    points_by_rib = {i: [] for i in range(len(rib_center_lines))}
    max_len = math.sqrt(wing_shape.BoundBox.XLength**2 + wing_shape.BoundBox.YLength**2) * 2
    z = z_min
    slice_count = 0
    while z <= z_max + 1e-6:
        wires = _get_wing_wires_at_z(wing_shape, z)
        if wires:
            for idx, line in enumerate(rib_center_lines):
                res = _rib_plane_intersect_z_plane(line, plane_normal, z)
                if res is None:
                    continue
                p0, d = res
                seg = _clip_line_to_wing_wires(p0, d, wires, max_len)
                if seg is None:
                    continue
                p1, p2 = seg
                mid = (p1 + p2) * 0.5
                points_by_rib[idx].append(mid)
        z += z_step
        slice_count += 1
    total = sum(len(v) for v in points_by_rib.values())
    print(f"Collected {total} midpoints over {slice_count} slices for {len(rib_center_lines)} ribs.")
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
