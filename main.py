# main.py
import FreeCAD
import FreeCADGui
import Part
from slab_utils import LWInfillParams, create_cut_result
from point_utils import collect_rib_midpoints, create_rib_wires
from prism_utils import create_bridges_trimmed_to_wing

# Force GUI if running from command line
if not FreeCAD.GuiUp:
    FreeCADGui.showMainWindow()
    import time
    time.sleep(1)  # allow GUI to initialise

if __name__ == "__main__":
    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)

    wing = doc.getObject("Pad")
    if not wing:
        raise RuntimeError("Object 'Pad' not found in document.")

    params = LWInfillParams(
        nozzle_diameter  = 0.4,
        rib_spacing      = 20.0,
        rib_width        = 0.1,
        rib_angle        = 30.0,
        grid_orientation = 0.0,
        primary_dir      = FreeCAD.Vector(0, 0, 1),
        construction_plane = 'XZ'
    )

    # 1. Create cut result (wing minus ribs)
    cut_result, all_center_lines, all_rib_solids = create_cut_result(wing.Shape, params, doc)
    print("Cut created.")

    # --- Visualise cut wing (semi‑transparent) ---
    cut_obj = doc.addObject("Part::Feature", "CutWing")
    # cut_obj.Shape = cut_result
    # doc.recompute()
    if FreeCAD.GuiUp:
        vp = FreeCADGui.ActiveDocument.getObject(cut_obj.Name)
        if vp:
            vp.Transparency = 80
            FreeCADGui.SendMsgToActiveView("ViewFit")
    else:
        print("No GUI – transparency not set.")

    # --- Visualise original rib centre lines (red) ---
    rib_lines_compound = Part.Compound([line for line in all_center_lines])
    # obj_lines = doc.addObject("Part::Feature", "RibCentreLines")
    # obj_lines.Shape = rib_lines_compound
    # if FreeCAD.GuiUp:
    #     line_view = FreeCADGui.ActiveDocument.getObject(obj_lines.Name)
    #     if line_view:
    #         line_view.LineColor = (1.0, 0.0, 0.0)
    doc.recompute()
    print(f"Visualized {len(all_center_lines)} rib centre lines.")

    # 2. Collect midpoints (handles holes → returns pieces)
    bb = wing.Shape.BoundBox
    z_min = bb.ZMin + 0.5
    z_max = bb.ZMax - 0.5
    z_step = 10.0
    print(f"Collecting midpoints from z={z_min:.1f} to {z_max:.1f}, step={z_step:.1f}...")
    pieces_by_rib = collect_rib_midpoints(
        wing_shape=wing.Shape,
        rib_center_lines=all_center_lines,
        plane_normal=params.plane_normal,
        z_min=z_min,
        z_max=z_max,
        z_step=z_step,
    )

    # 3. Create rib wires (green polylines) from the collected pieces
    create_rib_wires(pieces_by_rib, doc)

    # 4. Collect all points for a global point cloud
    all_points = []
    for pieces in pieces_by_rib.values():
        for piece_pts in pieces:
            all_points.extend(piece_pts)
    if all_points:
        vertices = [Part.Vertex(p) for p in all_points]
        compound = Part.Compound(vertices)
        obj = doc.addObject("Part::Feature", "AllPoints")
        # obj.Shape = compound
        # doc.recompute()
        print(f"Added 'AllPoints' with {len(all_points)} points.")

    # 5. Create trimmed bridges (one bridge per rib piece)
    bridges = create_bridges_trimmed_to_wing(
        data_by_rib=pieces_by_rib,
        rib_center_lines=all_center_lines,
        plane_normal=params.plane_normal,
        rib_width=params.rib_width,
        bridge_height=0.50,
        wing_shape=wing.Shape,
        extend_length=5.0,
        doc=doc
    )

    # 6. Merge cut wing with bridges to create final solid
    if bridges and hasattr(bridges, 'Shape') and not bridges.Shape.isNull():
        final_solid = cut_result.fuse(bridges.Shape)
        if final_solid and not final_solid.isNull():
            final_obj = doc.addObject("Part::Feature", "WingWithBridges")
            final_obj.Shape = final_solid
            doc.recompute()
            if FreeCAD.GuiUp:
                final_obj.ViewObject.Visibility = True
                final_obj.ViewObject.ShapeColor = (0.8, 0.8, 0.8)
                FreeCADGui.SendMsgToActiveView("ViewFit")
            print("Final solid (wing + bridges) created successfully.")
        else:
            print("Fusion resulted in null shape.")
    else:
        print("Bridges missing or invalid – final solid not created.")

    # Save the document so you can open it later
    doc.save()
    print("Document saved.")