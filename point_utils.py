# point_utils.py
import FreeCAD
import Part
import math
from collections import defaultdict

# ------------------------------------------------------------
# helpers (unchanged)
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
# wing cross‑section at Z (unchanged – for interior points)
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
# rib plane ∩ Z plane → line (for interior sampling)
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
# clip line to wing cross‑section → segment (for interior points)
# ------------------------------------------------------------
def _clip_line_to_wing_wires(point_on_line, line_dir, wires, max_len=2000.0):
    """
    Return a list of tuples (p1, p2, mid) for all inside segments
    where the line passes through the wing cross-section.
    """
    if not wires:
        return []
    p = point_on_line
    d = line_dir
    seg_start = FreeCAD.Vector(p.x - d.x*max_len, p.y - d.y*max_len, p.z)
    seg_end   = FreeCAD.Vector(p.x + d.x*max_len, p.y + d.y*max_len, p.z)
    try:
        test_edge = Part.makeLine(seg_start, seg_end)
    except Exception:
        return []

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
        return []
    t_vals.sort()
    segments = []
    for k in range(len(t_vals)-1):
        t1 = t_vals[k]
        t2 = t_vals[k+1]
        if abs(t2-t1) < 1e-6:
            continue
        t_mid = (t1+t2)/2.0
        mid_pt = FreeCAD.Vector(p.x + d.x*t_mid, p.y + d.y*t_mid, p.z)
        inside = False
        for wire in wires:
            try:
                face = Part.Face(wire)
                inside = face.isInside(mid_pt, 1e-3, True)
                if inside: break
            except:
                inside = wire.BoundBox.isInside(mid_pt)
                if inside: break
        if inside:
            p1 = FreeCAD.Vector(p.x + d.x*t1, p.y + d.y*t1, p.z)
            p2 = FreeCAD.Vector(p.x + d.x*t2, p.y + d.y*t2, p.z)
            mid = (p1 + p2) * 0.5
            segments.append((p1, p2, mid))
    return segments


# ------------------------------------------------------------
# MAIN COLLECTION (corrected edge‑case detection)
# ------------------------------------------------------------
def collect_rib_midpoints(wing_shape, rib_center_lines, plane_normal, z_min, z_max, z_step):
    from collections import defaultdict
    try:
        segments_per_rib = defaultdict(list)
        max_len = math.sqrt(wing_shape.BoundBox.XLength**2 + wing_shape.BoundBox.YLength**2) * 2

        # 1) Regular horizontal sampling
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
                    seg_list = _clip_line_to_wing_wires(p0, d, wires, max_len)
                    for i_seg, (p1, p2, mid) in enumerate(seg_list):
                        segments_per_rib[idx].append((z, p1, p2, mid, i_seg))
            z += z_step
            slice_count += 1

        # 2) Edge-case points from plane‑wing intersection
        edge_pts_per_rib = {}
        for idx, line in enumerate(rib_center_lines):
            start = line.Vertexes[0].Point
            end = line.Vertexes[-1].Point
            d = (end - start).normalize()
            n_rib = d.cross(plane_normal).normalize()
            plane = Part.Plane(start, n_rib)
            try:
                intersection = wing_shape.section(plane)
                if intersection and intersection.Edges:
                    vertices = []
                    for edge in intersection.Edges:
                        vertices.append(edge.Vertexes[0].Point)
                        vertices.append(edge.Vertexes[-1].Point)
                    if vertices:
                        min_z = min(p.z for p in vertices)
                        max_z = max(p.z for p in vertices)
                        z_tol = 1e-4
                        low_pts = [p for p in vertices if abs(p.z - min_z) <= z_tol]
                        high_pts = [p for p in vertices if abs(p.z - max_z) <= z_tol]
                        edge_pts = []
                        if low_pts:
                            c = FreeCAD.Vector(0,0,0)
                            for p in low_pts: c += p
                            edge_pts.append(c / len(low_pts))
                        if high_pts:
                            c = FreeCAD.Vector(0,0,0)
                            for p in high_pts: c += p
                            edge_pts.append(c / len(high_pts))
                        edge_pts_per_rib[idx] = edge_pts
            except Exception:
                pass

        # 3) Group segments into pieces and add edge-case points
        data_by_rib = {}
        for idx, entries in segments_per_rib.items():
            entries.sort(key=lambda x: x[0])   # by Z
            # group by segment index
            piece_dict = defaultdict(list)
            for z, p1, p2, mid, i_seg in entries:
                piece_dict[i_seg].append((z, mid))
            piece_data = []
            for piece_pts in piece_dict.values():
                pts_sorted = sorted(piece_pts, key=lambda x: x[0])
                midpoints = [p for (_, p) in pts_sorted]
                if midpoints:
                    zmin = min(p.z for p in midpoints)
                    zmax = max(p.z for p in midpoints)
                    piece_data.append((midpoints, zmin, zmax))
            # assign edge points to the correct piece
            if idx in edge_pts_per_rib:
                for ep in edge_pts_per_rib[idx]:
                    best_i = None
                    best_dist = float('inf')
                    for i, (_, zmin, zmax) in enumerate(piece_data):
                        if zmin - 1e-4 <= ep.z <= zmax + 1e-4:
                            best_i = i
                            break
                        dist = min(abs(ep.z - zmin), abs(ep.z - zmax))
                        if dist < best_dist:
                            best_dist = dist
                            best_i = i
                    if best_i is not None:
                        midpoints, _, _ = piece_data[best_i]
                        if not any(ep.isEqual(m, 1e-3) for m in midpoints):
                            midpoints.append(ep)
            # sort and deduplicate each piece
            pieces = []
            for midpoints, _, _ in piece_data:
                midpoints.sort(key=lambda p: p.z)
                uniq = []
                for p in midpoints:
                    if not uniq or not p.isEqual(uniq[-1], 1e-3):
                        uniq.append(p)
                pieces.append(uniq)
            if pieces:
                data_by_rib[idx] = pieces

        total_pieces = sum(len(v) for v in data_by_rib.values())
        print(f"Collected {total_pieces} rib pieces (including edge cases) over {slice_count} slices.")
        return data_by_rib if data_by_rib else {}
    except Exception as e:
        print(f"Fatal error in collect_rib_midpoints: {e}")
        import traceback
        traceback.print_exc()
        return {}
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


def create_rib_wires(data_by_rib, doc):
    wires = []
    for pieces in data_by_rib.values():
        for piece_pts in pieces:
            if len(piece_pts) < 2:
                continue
            pts_sorted = sorted(piece_pts, key=lambda p: p.z)
            edges = [Part.makeLine(pts_sorted[i], pts_sorted[i+1]) for i in range(len(pts_sorted)-1)]
            if edges:
                wire = Part.Wire(edges) if len(edges) > 1 else Part.Wire(edges[0])
                wires.append(wire)
    if wires:
        compound = Part.Compound(wires)
        obj = doc.addObject("Part::Feature", "RibWires")
        obj.Shape = compound
        if FreeCAD.GuiUp:
            obj.ViewObject.LineColor = (0.0, 1.0, 0.0)
            obj.ViewObject.LineWidth = 2
        doc.recompute()
        print(f"Created {len(wires)} rib wires.")
    else:
        print("No wires created.")