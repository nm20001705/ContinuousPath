# main.py
import FreeCAD
import FreeCADGui
import Part
from slab_utils import LWInfillParams, create_cut_result, generate_final_solid
from point_utils import collect_rib_midpoints, show_points_per_rib



if __name__ == "__main__":
    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)

    wing = doc.getObject("Pad")
    if not wing:
        raise RuntimeError("Object 'Pad' not found in document.")

    params = LWInfillParams(
        nozzle_diameter  = 0.4,
        rib_spacing      = 70.0,
        rib_width        = 0.1,
        rib_angle        = 30.0,
        grid_orientation = 0.0,
        primary_dir      = FreeCAD.Vector(0, 0, 1),
        construction_plane = 'XZ'
    )

    # 1. Create cut result
    cut_result, all_center_lines, all_rib_solids = create_cut_result(wing.Shape, params, doc)
    print("Cut created.")

    # --- Visualise cut wing ---
    cut_obj = doc.addObject("Part::Feature", "CutWing")
    cut_obj.Shape = cut_result
    doc.recompute()

    if FreeCAD.GuiUp:
        FreeCADGui.SendMsgToActiveView("ViewFit")
        view_provider = FreeCADGui.ActiveDocument.getObject(cut_obj.Name)
        if view_provider and hasattr(view_provider, "ViewObject"):
            view_provider.ViewObject.Transparency = 80
            print("Cut wing set to 80% transparency.")
        else:
            print("Could not set transparency – view object not ready.")
    else:
        print("No GUI – transparency not set.")

    # 2. Visualise rib centre lines
    rib_lines_compound = Part.Compound([line for line in all_center_lines])
    obj_lines = doc.addObject("Part::Feature", "RibCentreLines")
    obj_lines.Shape = rib_lines_compound
    if FreeCAD.GuiUp:
        line_view = FreeCADGui.ActiveDocument.getObject(obj_lines.Name)
        if line_view and hasattr(line_view, "ViewObject"):
            line_view.ViewObject.LineColor = (1.0, 0.0, 0.0)  # red
    doc.recompute()
    print(f"Visualized {len(all_center_lines)} rib centre lines.")

    # 3. Collect midpoints
    bb = wing.Shape.BoundBox
    z_min = bb.ZMin + 0.5
    z_max = bb.ZMax - 0.5
    z_step = 10
    print(f"Collecting midpoints from z={z_min:.1f} to {z_max:.1f}, step={z_step:.1f}...")
    points_by_rib = collect_rib_midpoints(
        wing_shape=wing.Shape,
        rib_center_lines=all_center_lines,
        plane_normal=params.plane_normal,
        z_min=bb.ZMin + 0.5,
        z_max=bb.ZMax - 0.5,
        z_step=z_step,
    )

    show_points_per_rib(points_by_rib, doc, mode='mid', prefix="RibMidpoints")
    # show_points_per_rib(points_by_rib, doc, mode='edge_cases', prefix="RibEdgeCases")

    # Combine all points for a global view
    all_points = []
    for data in points_by_rib.values():
        all_points.extend(data['mid'])
        all_points.extend(data['edge_cases'])
    if all_points:
        vertices = [Part.Vertex(p) for p in all_points]
        compound = Part.Compound(vertices)
        obj = doc.addObject("Part::Feature", "AllPoints")
        obj.Shape = compound
        doc.recompute()
        print(f"Added 'AllPoints' with {len(all_points)} points.")