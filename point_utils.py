# point_utils.py – as originally provided
import FreeCAD
import Part
import math
from collections import defaultdict

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------
def _cross(a, b):
    return FreeCAD.Vector(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x,
    )

def _dot(a, b):
    return a.x*b.x + a.y*b.y + a.z*b.z

def _norm(v):
    L = math.sqrt(v.x**2 + v.y**2 + v.z**2)
    if L < 1e-10:
        raise ValueError("Cannot normalize zero vector")
    return FreeCAD.Vector(v.x/L, v.y/L, v.z/L)

# ------------------------------------------------------------
# wing cross‑section at Z
# ------------------------------------------------------------
def _get_wing_wires_at_z(wing_shape, z):
    try:
        sliced = wing_shape.slice(FreeCAD.Vector(0, 0, 1), z)
        if not sliced:
            return []
        compounds = sliced if isinstance(sliced, list) else [sliced]
        wires = []
        for comp in compounds:
            if not hasattr(comp, 'Edges') or not comp.Edges:
                continue
            try:
                wires.append(Part.Wire(comp.Edges))
            except Exception:
                edges = list(comp.Edges)
                used = [False] * len(edges)
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
                            if e.Vertexes[0].Point.isEqual(cur_end, 1e-6):
                                loop.append(e)
                                used[j] = True
                                cur_end = e.Vertexes[-1].Point
                                changed = True
                                break
                            elif e.Vertexes[-1].Point.isEqual(cur_end, 1e-6):
                                rev = Part.Edge(e.Curve, e.Vertexes[-1].Point, e.Vertexes[0].Point)
                                loop.append(rev)
                                used[j] = True
                                cur_end = e.Vertexes[0].Point
                                changed = True
                                break
                    if len(loop) > 1:
                        try:
                            wires.append(Part.Wire(loop))
                        except Exception:
                            pass
        return wires
    except Exception:
        return []

# ------------------------------------------------------------
# rib plane ∩ Z plane → line
# ------------------------------------------------------------
def _rib_plane_intersect_z_plane(rib_center_line, plane_normal, z):
    start = rib_center_line.Vertexes[0].Point
    end   = rib_center_line.Vertexes[-1].Point
    try:
        d = _norm(FreeCAD.Vector(end.x-start.x, end.y-start.y, end.z-start.z))
    except ValueError:
        return None
    try:
        n_rib = _norm(_cross(d, plane_normal))
    except ValueError:
        return None
    Z = FreeCAD.Vector(0, 0, 1)
    line_dir_raw = _cross(n_rib, Z)
    try:
        line_dir = _norm(line_dir_raw)
    except ValueError:
        return None
    c = _dot(n_rib, start)
    rhs = c - n_rib.z * z
    if abs(n_rib.x) >= abs(n_rib.y):
        if abs(n_rib.x) < 1e-10:
            return None
        pt = FreeCAD.Vector(rhs / n_rib.x, 0.0, z)
    else:
        if abs(n_rib.y) < 1e-10:
            return None
        pt = FreeCAD.Vector(0.0, rhs / n_rib.y, z)
    return pt, line_dir

# ------------------------------------------------------------
# clip line to wing cross‑section → segment
# ------------------------------------------------------------
def _clip_line_to_wing_wires(point_on_line, line_dir, wires, max_len=2000.0):
    if not wires:
        return None
    p = point_on_line
    d = line_dir
    seg_start = FreeCAD.Vector(p.x - d.x*max_len, p.y - d.y*max_len, p.z)
    seg_end   = FreeCAD.Vector(p.x + d.x*max_len, p.y + d.y*max_len, p.z)
    try:
        test_edge = Part.makeLine(seg_start, seg_end)
    except Exception:
        return None
    t_vals = []
    for wire in wires:
        for edge in wire.Edges:
            try:
                dist, pts, _ = test_edge.distToShape(edge)
                if dist < 0.1:
                    for pt_pair in pts:
                        pt = pt_pair[0]
                        t = _dot(FreeCAD.Vector(pt.x-p.x, pt.y-p.y, 0.0), d)
                        t_vals.append(t)
            except Exception:
                continue
    if len(t_vals) < 2:
        return None
    t_vals.sort()
    best = None
    best_len = 0.0
    for k in range(len(t_vals)-1):
        t_mid = (t_vals[k] + t_vals[k+1]) / 2.0
        mid_pt = FreeCAD.Vector(p.x + d.x*t_mid, p.y + d.y*t_mid, p.z)
        inside = False
        for wire in wires:
            try:
                face = Part.Face(wire)
                inside = face.isInside(mid_pt, 1e-3, True)
            except Exception:
                inside = wire.BoundBox.isInside(mid_pt)
            if inside:
                break
        if inside:
            seg_len = t_vals[k+1] - t_vals[k]
            if seg_len > best_len:
                best_len = seg_len
                p1 = FreeCAD.Vector(p.x + d.x*t_vals[k], p.y + d.y*t_vals[k], p.z)
                p2 = FreeCAD.Vector(p.x + d.x*t_vals[k+1], p.y + d.y*t_vals[k+1], p.z)
                mid = FreeCAD.Vector((p1.x+p2.x)/2, (p1.y+p2.y)/2, p.z)
                best = (p1, p2, mid)
    return best

# ------------------------------------------------------------
# MAIN COLLECTION
# ------------------------------------------------------------
def collect_rib_midpoints(wing_shape, rib_center_lines, plane_normal, z_min, z_max, z_step, doc=None, vis=False):
    data_by_rib = defaultdict(lambda: {'mid': [], 'edge_cases': []})
    max_len = math.sqrt(wing_shape.BoundBox.XLength**2 + wing_shape.BoundBox.YLength**2) * 2

    # 1) Regular horizontal sampling – collect midpoints
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
                _, _, mid = seg
                data_by_rib[idx]['mid'].append(mid)
        z += z_step
        slice_count += 1

    # 2) Edge‑case points from plane‑wing intersection
    for idx, line in enumerate(rib_center_lines):
        start = line.Vertexes[0].Point
        end   = line.Vertexes[-1].Point
        d = (end - start).normalize()
        n_rib = d.cross(plane_normal).normalize()
        plane = Part.Plane(start, n_rib)

        try:
            intersection = wing_shape.section(plane)
            if not intersection or not intersection.Edges:
                continue

            vertices = []
            for edge in intersection.Edges:
                vertices.append(edge.Vertexes[0].Point)
                vertices.append(edge.Vertexes[-1].Point)

            if not vertices:
                continue

            min_z = min(p.z for p in vertices)
            max_z = max(p.z for p in vertices)
            z_tol = 1e-4
            low_pts  = [p for p in vertices if abs(p.z - min_z) <= z_tol]
            high_pts = [p for p in vertices if abs(p.z - max_z) <= z_tol]

            unique_low = []
            for p in low_pts:
                if not any(p == q for q in unique_low):
                    unique_low.append(p)
            unique_high = []
            for p in high_pts:
                if not any(p == q for q in unique_high):
                    unique_high.append(p)

            data_by_rib[idx]['edge_cases'] = unique_low + unique_high

        except Exception:
            # Fallback: use min/max from midpoints
            if data_by_rib[idx]['mid']:
                mid_pts = data_by_rib[idx]['mid']
                min_z = min(p.z for p in mid_pts)
                max_z = max(p.z for p in mid_pts)
                z_tol = 1e-4
                low_pts  = [p for p in mid_pts if abs(p.z - min_z) <= z_tol]
                high_pts = [p for p in mid_pts if abs(p.z - max_z) <= z_tol]
                unique_low = []
                for p in low_pts:
                    if not any(p == q for q in unique_low):
                        unique_low.append(p)
                unique_high = []
                for p in high_pts:
                    if not any(p == q for q in unique_high):
                        unique_high.append(p)
                data_by_rib[idx]['edge_cases'] = unique_low + unique_high

    # 3) Merge edge cases into midpoints
    for idx, data in data_by_rib.items():
        if not data['edge_cases']:
            continue
        edge_pts = data['edge_cases']
        z_groups = {}
        for p in edge_pts:
            z_key = round(p.z, 4)
            z_groups.setdefault(z_key, []).append(p)
        for z_key, pts in z_groups.items():
            if len(pts) == 1:
                centroid = pts[0]
            else:
                centroid = FreeCAD.Vector(0,0,0)
                for p in pts:
                    centroid += p
                centroid /= len(pts)
            already = any(centroid.isEqual(m, 1e-3) for m in data['mid'])
            if not already:
                data['mid'].append(centroid)

    total_mid = sum(len(v['mid']) for v in data_by_rib.values())
    total_edge = sum(len(v['edge_cases']) for v in data_by_rib.values())
    print(f"Collected {total_mid} midpoints over {slice_count} slices.")
    if vis:
        all_points = []
        for data in data_by_rib.values():
            all_points.extend(data['mid'])
            all_points.extend(data['edge_cases'])
        if all_points:
            from viz_utils import show_midpoints
            show_midpoints(all_points, doc, point_size=5)

    return data_by_rib

# ------------------------------------------------------------
# VISUALISATION
# ------------------------------------------------------------
def show_points_per_rib(data_by_rib, doc, mode='mid', prefix='RibPoints'):
    count = 0
    for idx, data in data_by_rib.items():
        if mode == 'mid':
            pts = data['mid']
        elif mode == 'edge_cases':
            pts = data['edge_cases']
        elif mode == 'all':
            pts = data['mid'] + data['edge_cases']
        else:
            raise ValueError("mode must be 'mid', 'edge_cases', or 'all'")
        if not pts:
            continue
        vertices = [Part.Vertex(p) for p in pts]
        compound = Part.Compound(vertices)
        obj = doc.addObject("Part::Feature", f"{prefix}_{idx}")
        obj.Shape = compound
        if FreeCAD.GuiUp:
            obj.ViewObject.PointSize = 5
            obj.ViewObject.DisplayMode = "Points"
        count += 1
    doc.recompute()
    total = sum(len(data['mid']) if mode=='mid' else len(data['edge_cases']) if mode=='edge_cases' else len(data['mid'])+len(data['edge_cases']) for data in data_by_rib.values())
    print(f"Visualized {count} ribs with {total} points (mode={mode}).")

def create_rib_wires(data_by_rib, doc, vis=False):
    wires = []
    for idx, data in data_by_rib.items():
        pts = data['mid']
        if len(pts) < 2:
            continue
        pts_sorted = sorted(pts, key=lambda p: p.z)
        edges = [Part.makeLine(pts_sorted[i], pts_sorted[i+1]) for i in range(len(pts_sorted)-1)]
        if edges:
            wire = Part.Wire(edges) if len(edges) > 1 else Part.Wire(edges[0])
            wires.append(wire)
    if vis and wires:
        from viz_utils import show_rib_wires
        show_rib_wires(wires, doc)
    return wires