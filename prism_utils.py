# prism_utils.py
import FreeCAD
import Part

def create_rectangular_prism_bridge(data_by_rib, rib_center_lines, plane_normal, rib_width, bridge_height, doc=None):
    """
    Create a solid rectangular prism around each rib wire by offsetting the
    centreline in two perpendicular directions: rib normal and plane normal.

    Parameters:
        data_by_rib      : output of collect_rib_midpoints (contains 'mid' points)
        rib_center_lines : original rib centre lines (to compute rib normals)
        plane_normal     : construction plane normal (extrusion direction)
        rib_width        : thickness of the rib cut (gap width)
        bridge_height    : height of the bridge perpendicular to the rib plane
        doc              : FreeCAD document
    """
    if doc is None:
        doc = FreeCAD.ActiveDocument

    # Pre‑compute rib normals from centre lines
    rib_normals = []
    for line in rib_center_lines:
        start = line.Vertexes[0].Point
        end   = line.Vertexes[-1].Point
        rib_dir = (end - start).normalize()
        n_rib = rib_dir.cross(plane_normal).normalize()
        rib_normals.append(n_rib)

    all_bridges = []

    for idx, data in data_by_rib.items():
        pts = data['mid']
        if len(pts) < 2:
            continue
        pts_sorted = sorted(pts, key=lambda p: p.z)   # bottom → top
        n_rib = rib_normals[idx]
        ext_dir = plane_normal.normalize()
        half_w = rib_width / 2.0
        half_h = bridge_height / 2.0

        # For each segment, build a prism by connecting the offset quadrilaterals
        seg_prisms = []
        for i in range(len(pts_sorted)-1):
            p1 = pts_sorted[i]
            p2 = pts_sorted[i+1]

            # Eight vertices (four at p1, four at p2)
            offs = [
                ( half_w,  half_h),
                ( half_w, -half_h),
                (-half_w, -half_h),
                (-half_w,  half_h),
            ]
            v1 = [p1 + n_rib * w + ext_dir * h for (w, h) in offs]
            v2 = [p2 + n_rib * w + ext_dir * h for (w, h) in offs]

            # Build six faces (4 side faces + 2 end caps)
            faces = []
            # Side faces (quadrilaterals)
            for j in range(4):
                k = (j+1) % 4
                face_vertices = [v1[j], v2[j], v2[k], v1[k]]
                wire = Part.makePolygon(face_vertices + [face_vertices[0]])
                faces.append(Part.Face(wire))
            # End faces at p1 and p2
            wire1 = Part.makePolygon(v1 + [v1[0]])
            faces.append(Part.Face(wire1))
            wire2 = Part.makePolygon(v2 + [v2[0]])
            faces.append(Part.Face(wire2))

            # Solid from faces
            try:
                shell = Part.makeShell(faces)
                if not shell.isNull():
                    solid = Part.makeSolid(shell)
                    if not solid.isNull():
                        seg_prisms.append(solid)
            except Exception as e:
                print(f"Error creating prism for rib {idx} segment {i}: {e}")

        if seg_prisms:
            # Fuse all prisms of this rib into one solid
            fused = seg_prisms[0]
            for pr in seg_prisms[1:]:
                fused = fused.fuse(pr)
            all_bridges.append(fused)

    if all_bridges:
        compound = Part.Compound(all_bridges)
        obj = doc.addObject("Part::Feature", "Bridges")
        obj.Shape = compound
        if FreeCAD.GuiUp:
            obj.ViewObject.ShapeColor = (0.0, 0.8, 0.0)
            obj.ViewObject.Transparency = 30
        doc.recompute()
        print(f"Created {len(all_bridges)} bridges (one per rib).")
        return obj
    else:
        print("No bridges created.")
        return None