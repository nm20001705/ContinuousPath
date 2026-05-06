import FreeCAD
import Part
import math
from collections import defaultdict


# ============================================================
# HELPERS
# ============================================================

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


# ============================================================
# WING CROSS-SECTION AT Z
# ============================================================

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
                    if len(loop) > 1:
                        try:
                            wires.append(Part.Wire(loop))
                        except Exception:
                            pass
        return wires
    except Exception:
        return []


# ============================================================
# RIB FACE PLANE ∩ Z PLANE → LINE
# ============================================================

def _face_plane_intersect_z_plane(ref_point, rib_dir, plane_normal, z):
    """
    Intersect the plane defined by (ref_point, rib_dir, plane_normal)
    with z = const.

    ref_point  : a point on the rib face plane (offset from center by ±rib_width/2)
    rib_dir    : direction along the rib (unit vector)
    plane_normal: construction plane extrusion direction

    Rib face plane normal: n_face = rib_dir × plane_normal
    Intersection line direction: line_dir = n_face × Z

    Returns (point_on_line, line_dir) or None.
    """
    try:
        d = _norm(rib_dir)
    except ValueError:
        return None
    try:
        n_face = _norm(_cross(d, plane_normal))
    except ValueError:
        return None

    Z = FreeCAD.Vector(0, 0, 1)
    try:
        line_dir = _norm(_cross(n_face, Z))
    except ValueError:
        return None  # face plane parallel to z = const

    # Find point on intersection line:
    # n_face · p = n_face · ref_point, with p.z = z
    c   = _dot(n_face, ref_point)
    rhs = c - n_face.z * z

    if abs(n_face.x) >= abs(n_face.y):
        if abs(n_face.x) < 1e-10:
            return None
        pt = FreeCAD.Vector(rhs / n_face.x, 0.0, z)
    else:
        if abs(n_face.y) < 1e-10:
            return None
        pt = FreeCAD.Vector(0.0, rhs / n_face.y, z)

    return pt, line_dir


# ============================================================
# CLIP LINE TO WING CROSS-SECTION → SEGMENT
# ============================================================

def _clip_line_to_wing_wires(point_on_line, line_dir, wires, max_len=2000.0):
    """
    Clip infinite line to wing cross-section wires.
    Returns (p1, p2, mid) — exact boundary endpoints and midpoint,
    or None if no valid inside segment found.
    """
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

    for k in range(len(t_vals) - 1):
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
                p1  = FreeCAD.Vector(p.x + d.x*t_vals[k],   p.y + d.y*t_vals[k],   p.z)
                p2  = FreeCAD.Vector(p.x + d.x*t_vals[k+1], p.y + d.y*t_vals[k+1], p.z)
                mid = FreeCAD.Vector((p1.x+p2.x)/2, (p1.y+p2.y)/2, p.z)
                best = (p1, p2, mid)

    return best


# ============================================================
# RIB NORMAL COMPUTATION
# ============================================================

def _rib_normal(rib_center_line, plane_normal):
    """
    Compute the rib face normal: n = rib_dir × plane_normal.
    This points perpendicular to the rib face, in the slab thickness direction.
    """
    start = rib_center_line.Vertexes[0].Point
    end   = rib_center_line.Vertexes[-1].Point
    try:
        d = _norm(FreeCAD.Vector(end.x-start.x, end.y-start.y, end.z-start.z))
        return _norm(_cross(d, plane_normal))
    except ValueError:
        return None


# ============================================================
# MAIN COLLECTION
# ============================================================

def collect_rib_midpoints(wing_shape, rib_center_lines, plane_normal,
                          z_min, z_max, z_step, rib_width):
    """
    For each rib and each Z slice, compute midpoints on the two large
    face planes of the slab (offset ±rib_width/2 from center along
    the rib normal). Each face gives its own segment and midpoint.

    Returns
    -------
    dict {rib_index: {'face_a': [FreeCAD.Vector, ...],
                      'face_b': [FreeCAD.Vector, ...]}}
    where face_a and face_b are the two slab face midpoint sequences,
    sorted by Z, ready for surface construction.
    """
    bb = wing_shape.BoundBox
    max_len = math.sqrt(bb.XLength**2 + bb.YLength**2) * 2
    half_w = rib_width / 2.0

    data_by_rib = {}
    for idx in range(len(rib_center_lines)):
        data_by_rib[idx] = {'face_a': [], 'face_b': []}

    # Pre-compute rib direction and normal for each rib
    rib_info = []
    for line in rib_center_lines:
        start = line.Vertexes[0].Point
        end   = line.Vertexes[-1].Point
        try:
            d = _norm(FreeCAD.Vector(end.x-start.x, end.y-start.y, end.z-start.z))
            n = _norm(_cross(d, plane_normal))
        except ValueError:
            d, n = None, None
        rib_info.append({'start': start, 'dir': d, 'normal': n})

    # Slice loop
    z = z_min
    slice_count = 0
    while z <= z_max + 1e-6:
        wires = _get_wing_wires_at_z(wing_shape, z)
        if wires:
            for idx, info in enumerate(rib_info):
                if info['dir'] is None or info['normal'] is None:
                    continue

                d = info['dir']
                n = info['normal']
                start = info['start']

                # Move the reference point to the current Z level along the rib.
                # The rib plane is infinite so any point on it works — but using
                # a point near the wing cross-section avoids numerical issues.
                # Project start onto this Z by stepping along rib direction.
                if abs(d.z) > 1e-6:
                    t_z = (z - start.z) / d.z
                    center_at_z = FreeCAD.Vector(
                        start.x + d.x * t_z,
                        start.y + d.y * t_z,
                        z,
                    )
                else:
                    center_at_z = FreeCAD.Vector(start.x, start.y, z)

                # Two reference points — one on each face at this Z
                ref_a = FreeCAD.Vector(
                    center_at_z.x + n.x * half_w,
                    center_at_z.y + n.y * half_w,
                    center_at_z.z + n.z * half_w,
                )
                ref_b = FreeCAD.Vector(
                    center_at_z.x - n.x * half_w,
                    center_at_z.y - n.y * half_w,
                    center_at_z.z - n.z * half_w,
                )

                for ref, face_key in ((ref_a, 'face_a'), (ref_b, 'face_b')):
                    res = _face_plane_intersect_z_plane(ref, d, plane_normal, z)
                    if res is None:
                        continue
                    pt, line_dir = res
                    seg = _clip_line_to_wing_wires(pt, line_dir, wires, max_len)
                    if seg is None:
                        continue
                    _, _, mid = seg
                    data_by_rib[idx][face_key].append(mid)

        z += z_step
        slice_count += 1

    # Sort each face sequence by Z
    for idx in data_by_rib:
        for face_key in ('face_a', 'face_b'):
            data_by_rib[idx][face_key].sort(key=lambda p: p.z)

    total = sum(
        len(d['face_a']) + len(d['face_b'])
        for d in data_by_rib.values()
    )
    print(f"Collected {total} face midpoints over {slice_count} slices "
          f"for {len(rib_center_lines)} ribs.")
    return data_by_rib


# ============================================================
# SURFACE CONSTRUCTION
# ============================================================

def create_rib_surfaces(data_by_rib, doc):
    """
    For each rib, build a ruled surface between face_a and face_b
    point sequences (both sorted bottom→top).
    Uses Part.makeRuledSurface which handles non-planar wires.
    """
    surfaces = []

    for idx, data in data_by_rib.items():
        pts_a = data['face_a']  # bottom→top
        pts_b = data['face_b']  # bottom→top

        if len(pts_a) < 2 or len(pts_b) < 2:
            print(f"  [skip] rib {idx}: not enough points "
                  f"(face_a={len(pts_a)}, face_b={len(pts_b)})")
            continue

        try:
            edges_a = [Part.makeLine(pts_a[i], pts_a[i+1])
                       for i in range(len(pts_a)-1)]
            wire_a = Part.Wire(edges_a)

            edges_b = [Part.makeLine(pts_b[i], pts_b[i+1])
                       for i in range(len(pts_b)-1)]
            wire_b = Part.Wire(edges_b)

            surface = Part.makeRuledSurface(wire_a, wire_b)
            if surface and not surface.isNull():
                surfaces.append(surface)
            else:
                print(f"  [warn] rib {idx}: ruled surface is null")
        except Exception as e:
            print(f"  [warn] rib {idx}: surface failed: {e}")

    if surfaces:
        compound = Part.Compound(surfaces)
        obj = doc.addObject("Part::Feature", "RibSurfaces")
        obj.Shape = compound
        doc.recompute()
        print(f"Created {len(surfaces)} rib surfaces.")
    else:
        print("No rib surfaces created.")
    return surfaces


# ============================================================
# VISUALISATION
# ============================================================

def show_points_per_rib(data_by_rib, doc, mode='both', prefix='RibPoints'):
    """
    mode: 'face_a'  → face A midpoints only
          'face_b'  → face B midpoints only
          'both'    → both faces combined
    """
    count = 0
    for idx, data in data_by_rib.items():
        if mode == 'face_a':
            pts = data['face_a']
        elif mode == 'face_b':
            pts = data['face_b']
        else:
            pts = data['face_a'] + data['face_b']
        if not pts:
            continue
        compound = Part.Compound([Part.Vertex(p) for p in pts])
        obj = doc.addObject("Part::Feature", f"{prefix}_{idx}")
        obj.Shape = compound
        count += 1
    doc.recompute()
    total = sum(
        len(d['face_a']) + len(d['face_b']) for d in data_by_rib.values()
    )
    print(f"Visualized {count} ribs ({total} points, mode='{mode}').")