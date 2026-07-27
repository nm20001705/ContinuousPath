# prism_utils.py - Mesh-optimized version with correct method names

import FreeCAD
import Part
import Mesh
import MeshPart

def extend_rib_wire(pts_sorted, extend_len):
    if len(pts_sorted) < 2:
        return pts_sorted
    dir0 = (pts_sorted[1] - pts_sorted[0]).normalize()
    p_before = pts_sorted[0] - dir0 * extend_len
    dir_last = (pts_sorted[-1] - pts_sorted[-2]).normalize()
    p_after = pts_sorted[-1] + dir_last * extend_len
    return [p_before] + pts_sorted + [p_after]

def create_bridges_trimmed_to_wing_mesh(data_by_rib, rib_center_lines, plane_normal,
                                        rib_width, bridge_height, wing_mesh,
                                        extend_length=5.0, doc=None, vis=False):
    if doc is None:
        doc = FreeCAD.ActiveDocument

    rib_normals = []
    for line in rib_center_lines:
        start = line.Vertexes[0].Point
        end   = line.Vertexes[-1].Point
        rib_dir = (end - start).normalize()
        n_rib = rib_dir.cross(plane_normal).normalize()
        rib_normals.append(n_rib)

    all_bridge_meshes = []

    for idx, data in data_by_rib.items():
        pts = data['mid']
        if len(pts) < 2:
            continue
        pts_sorted = sorted(pts, key=lambda p: p.z)
        pts_extended = extend_rib_wire(pts_sorted, extend_length)

        n_rib = rib_normals[idx]
        ext_dir = plane_normal.normalize()

        seg_meshes = []
        for i in range(len(pts_extended)-1):
            p1 = pts_extended[i]
            p2 = pts_extended[i+1]

            mid = (p1 + p2) * 0.5
            length = (p2 - p1).Length
            if length < 1e-6:
                continue
            local_x = (p2 - p1).normalize()
            local_y = n_rib
            local_z = ext_dir

            box = Part.makeBox(length, rib_width, bridge_height)
            mat = FreeCAD.Matrix()
            mat.A11, mat.A12, mat.A13 = local_x.x, local_y.x, local_z.x
            mat.A21, mat.A22, mat.A23 = local_x.y, local_y.y, local_z.y
            mat.A31, mat.A32, mat.A33 = local_x.z, local_y.z, local_z.z
            mat.A44 = 1.0
            box.Placement = FreeCAD.Placement(mid, FreeCAD.Rotation(mat))
            box.translate(FreeCAD.Vector(-length/2, -rib_width/2, -bridge_height/2))

            mesh = MeshPart.meshFromShape(box, LinearDeflection=0.5, AngularDeflection=0.5)
            seg_meshes.append(mesh)

        if seg_meshes:
            # Use "unite" instead of "union"
            rib_bridge = seg_meshes[0]
            for sm in seg_meshes[1:]:
                rib_bridge = rib_bridge.unite(sm)      # <-- changed
                rib_bridge.clean()
                rib_bridge.removeDuplicatedPoints()
                rib_bridge.removeDuplicatedFacets()
            trimmed = rib_bridge.intersect(wing_mesh)
            if trimmed and trimmed.CountFacets > 0:
                all_bridge_meshes.append(trimmed)

    if not all_bridge_meshes:
        return None

    final_bridges = all_bridge_meshes[0]
    for bm in all_bridge_meshes[1:]:
        final_bridges = final_bridges.unite(bm)        # <-- changed
        final_bridges.clean()
        final_bridges.removeDuplicatedPoints()
        final_bridges.removeDuplicatedFacets()

    if vis:
        from viz_utils import show_mesh
        show_mesh(final_bridges, doc, "BridgesMesh", color=(0.0,0.8,0.0), transparency=30)

    return final_bridges