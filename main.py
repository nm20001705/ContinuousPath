# main.py
import FreeCAD
import FreeCADGui
import Part
from slab_utils import LWInfillParams, create_cut_result, merge_and_show_final
from point_utils import collect_rib_midpoints, create_rib_wires
from prism_utils import create_bridges_trimmed_to_wing
from viz_utils import show_final_solid, fit_view

if __name__ == "__main__":
    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)

    wing = doc.getObject("Pad")
    if not wing:
        raise RuntimeError("Object 'Pad' not found.")

    params = LWInfillParams(
        nozzle_diameter=0.4,
        rib_spacing=70.0,
        rib_width=0.1,
        rib_angle=30.0,
        grid_orientation=0.0,
        z_step = 10,
        primary_dir=FreeCAD.Vector(0, 0, 1),
        construction_plane='XZ',
        vis_cut_wing=False,
        vis_centre_lines=False,
        vis_midpoints=False,
        vis_wires=False,
        vis_bridges=False,
        vis_final_solid=True
    )

    # 1. Create cut result (cut wing is shown internally if vis_cut_wing is True)
    cut_result, all_center_lines, _ = create_cut_result(wing.Shape, params, doc, vis=params.vis_cut_wing)
    print("Cut created.")

    # 2. Visualise rib centre lines (if flag is True)
    if params.vis_centre_lines:
        from viz_utils import show_rib_centre_lines
        show_rib_centre_lines(all_center_lines, doc)

    # 3. Collect midpoints and optionally visualise them
    bb = wing.Shape.BoundBox
    z_min = bb.ZMin + 0.5
    z_max = bb.ZMax - 0.5
    points_by_rib = collect_rib_midpoints(
        wing_shape=wing.Shape,
        rib_center_lines=all_center_lines,
        plane_normal=params.plane_normal,
        z_min=z_min, z_max=z_max, z_step=params.z_step,
        doc=doc, vis=params.vis_midpoints
    )

    # 4. Create rib wires and optionally visualise them
    wires = create_rib_wires(points_by_rib, doc, vis=params.vis_wires)

    # 5. Create bridges and optionally visualise them
    bridges = create_bridges_trimmed_to_wing(
        data_by_rib=points_by_rib,
        rib_center_lines=all_center_lines,
        plane_normal=params.plane_normal,
        rib_width=params.rib_width,
        bridge_height=0.5,
        wing_shape=wing.Shape,
        extend_length=5.0,
        doc=doc,
        vis=params.vis_bridges
    )

    # 6. Merge and show final solid
    if bridges:
        final_solid = merge_and_show_final(cut_result, bridges, doc, vis=params.vis_final_solid)
        print("Final solid (wing + bridges) created successfully.")
    else:
        print("No bridges created – final solid not produced.")

    fit_view(doc)
    doc.save()
    print("Document saved.")