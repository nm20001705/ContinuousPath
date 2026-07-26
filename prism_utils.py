# prism_utils.py
import FreeCAD
import Part

from perf_utils import timed

def extend_rib_wire(pts_sorted, extend_len):
    """Extend the first and last segment outward, return new list with two extra points."""
    if len(pts_sorted) < 2:
        return pts_sorted
    dir0 = (pts_sorted[1] - pts_sorted[0]).normalize()
    p_before = pts_sorted[0] - dir0 * extend_len
    dir_last = (pts_sorted[-1] - pts_sorted[-2]).normalize()
    p_after = pts_sorted[-1] + dir_last * extend_len
    return [p_before] + pts_sorted + [p_after]

@timed
def create_bridges_trimmed_to_wing(data_by_rib, rib_center_lines, plane_normal,
                                   rib_width, bridge_height, wing_shape,
                                   extend_length=5.0, doc=None, vis=False):
    """
    For each rib, create a rectangular prism along the rib wire, extend it,
    then intersect with the wing solid. Returns a compound of bridges.
    """
    if doc is None:
        doc = FreeCAD.ActiveDocument

    # Quick exit if no data
    if not data_by_rib or not any(v['mid'] for v in data_by_rib.values()):
        print("No bridge points – skipping bridge generation.")
        return None

    # Pre‑compute rib normals from centre lines
    rib_normals = []
    for line in rib_center_lines:
        start = line.Vertexes[0].Point
        end   = line.Vertexes[-1].Point
        rib_dir = (end - start).normalize()
        n_rib = rib_dir.cross(plane_normal).normalize()
        rib_normals.append(n_rib)

    all_bridges = []  # collect all trimmed bridge solids
    total_segments = 0
    total_prisms = 0
    ribs_with_bridges = 0

    for idx, data in data_by_rib.items():
        pts = data['mid']
        if len(pts) < 2:
            continue
        pts_sorted = sorted(pts, key=lambda p: p.z)
        pts_extended = extend_rib_wire(pts_sorted, extend_length)
        total_segments += len(pts_extended) - 1

        n_rib = rib_normals[idx]
        ext_dir = plane_normal.normalize()
        half_w = rib_width / 2.0
        half_h = bridge_height / 2.0

        seg_prisms = []
        for i in range(len(pts_extended)-1):
            p1 = pts_extended[i]
            p2 = pts_extended[i+1]
            offs = [
                ( half_w,  half_h),
                ( half_w, -half_h),
                (-half_w, -half_h),
                (-half_w,  half_h),
            ]
            v1 = [p1 + n_rib * w + ext_dir * h for (w, h) in offs]
            v2 = [p2 + n_rib * w + ext_dir * h for (w, h) in offs]

            faces = []
            for j in range(4):
                k = (j+1) % 4
                verts = [v1[j], v2[j], v2[k], v1[k]]
                wire = Part.makePolygon(verts + [verts[0]])
                faces.append(Part.Face(wire))
            wire1 = Part.makePolygon(v1 + [v1[0]])
            faces.append(Part.Face(wire1))
            wire2 = Part.makePolygon(v2 + [v2[0]])
            faces.append(Part.Face(wire2))

            try:
                shell = Part.makeShell(faces)
                if not shell.isNull():
                    solid = Part.makeSolid(shell)
                    if not solid.isNull():
                        seg_prisms.append(solid)
                        total_prisms += 1
            except Exception:
                continue

        if seg_prisms:
            # Fuse all segments of this rib into one solid
            fused = seg_prisms[0]
            for pr in seg_prisms[1:]:
                try:
                    fused = fused.fuse(pr)
                except:
                    pass

            # Trim the extended bridge to the wing shape
            try:
                trimmed = fused.common(wing_shape)
                if not trimmed.isNull() and trimmed.Volume > 1e-6:
                    all_bridges.append(trimmed)
                    ribs_with_bridges += 1
            except Exception:
                continue


    if all_bridges:
        compound = Part.Compound(all_bridges)
        if vis:
            from viz_utils import show_bridges
            show_bridges(compound, doc)
        return compound
    else:
        print("No bridges created.")
        return None
