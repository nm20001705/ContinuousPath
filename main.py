import FreeCAD
import FreeCADGui
import Part
from slab_utils import create_cut_result, merge_and_show_final, add_ellipse_holes_to_faces
from point_utils import collect_rib_midpoints, create_rib_wires
from prism_utils import create_bridges_trimmed_to_wing
from viz_utils import fit_view

PLANE_DEFS = {
    'XY': {'normal': FreeCAD.Vector(0, 0, 1),
           'axis_u': FreeCAD.Vector(1, 0, 0),
           'axis_v': FreeCAD.Vector(0, 1, 0)},
    'XZ': {'normal': FreeCAD.Vector(0, 1, 0),
           'axis_u': FreeCAD.Vector(1, 0, 0),
           'axis_v': FreeCAD.Vector(0, 0, 1)},
    'YZ': {'normal': FreeCAD.Vector(1, 0, 0),
           'axis_u': FreeCAD.Vector(0, 1, 0),
           'axis_v': FreeCAD.Vector(0, 0, 1)},
}

class LWInfillParams:
    def __init__(self, nozzle_diameter=0.4, wall_thickness=None,
                 rib_spacing=5.0, rib_width=None, rib_angle=30.0,
                 grid_orientation=0.0, primary_dir=None, z_step=10,
                 construction_plane='XZ',
                 vis_cut_wing=True,
                 vis_rib_centre_surfaces=True,
                 vis_centre_lines=True,
                 vis_midpoints=True,
                 vis_wires=True,
                 vis_bridges=True,
                 vis_rib_surface_segments=True,
                 vis_final_solid=True):
        if construction_plane not in PLANE_DEFS:
            raise ValueError(f"construction_plane must be one of {list(PLANE_DEFS.keys())}")
        self.nozzle_diameter = nozzle_diameter
        self.wall_thickness = wall_thickness or nozzle_diameter
        self.rib_spacing = rib_spacing
        self.rib_width = rib_width or (2.0 * nozzle_diameter)
        self.rib_angle = rib_angle
        self.grid_orientation = grid_orientation
        self.construction_plane = construction_plane
        pdef = PLANE_DEFS[construction_plane]
        self.plane_normal = pdef['normal']
        self.plane_axis_u = pdef['axis_u']
        self.plane_axis_v = pdef['axis_v']
        self.primary_dir = self._project_primary(primary_dir)
        self.z_step = z_step
        self.vis_cut_wing = vis_cut_wing
        self.vis_centre_lines = vis_centre_lines
        self.vis_midpoints = vis_midpoints
        self.vis_wires = vis_wires
        self.vis_bridges = vis_bridges
        self.vis_rib_surface_segments = vis_rib_surface_segments
        self.vis_rib_centre_surfaces = vis_rib_centre_surfaces
        self.vis_final_solid = vis_final_solid

    def _project_primary(self, pd):
        if pd is None:
            return None
        n = self.plane_normal
        dot = pd.x*n.x + pd.y*n.y + pd.z*n.z
        proj = FreeCAD.Vector(pd.x - dot*n.x, pd.y - dot*n.y, pd.z - dot*n.z)
        L = proj.Length
        if L < 1e-6:
            print(f"[warn] primary_dir parallel to {self.construction_plane} normal – auto-detect.")
            return None
        return proj.normalize()

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
        z_step=10,
        primary_dir=FreeCAD.Vector(0, 0, 1),
        construction_plane='XZ',
        vis_cut_wing=False,
        vis_centre_lines=False,
        vis_midpoints=False,
        vis_wires=False,
        vis_bridges=False,
        vis_rib_centre_surfaces=True,
        vis_rib_surface_segments=True,
        vis_final_solid=True,
    )

    # 1. Create cut result
    cut_result, all_center_lines, valid_ribs = create_cut_result(
        wing.Shape, params, doc, vis=params.vis_cut_wing
    )
    print("Cut created.")

    # 2. Rib centre lines
    if params.vis_centre_lines:
        from viz_utils import show_rib_centre_lines
        show_rib_centre_lines(all_center_lines, doc)

    # 3. Collect midpoints
    bb = wing.Shape.BoundBox
    points_by_rib = collect_rib_midpoints(
        wing_shape=wing.Shape,
        rib_center_lines=all_center_lines,
        plane_normal=params.plane_normal,
        z_min=bb.ZMin + 0.5,
        z_max=bb.ZMax - 0.5,
        z_step=params.z_step,
        doc=doc,
        vis=params.vis_midpoints,
    )

    # 4. Rib wires
    wires = create_rib_wires(points_by_rib, doc, vis=params.vis_wires)

    # 5. Bridges
    bridges = create_bridges_trimmed_to_wing(
        data_by_rib=points_by_rib,
        rib_center_lines=all_center_lines,
        plane_normal=params.plane_normal,
        rib_width=params.rib_width,
        bridge_height=0.5,
        wing_shape=wing.Shape,
        extend_length=5.0,
        doc=doc,
        vis=params.vis_bridges,
    )

    # 6. Final solid
    if bridges:
        final_solid = merge_and_show_final(cut_result, bridges, doc, vis=params.vis_final_solid)
        print("Final solid (wing + bridges) created successfully.")
    else:
        print("No bridges created – final solid not produced.")

    # 7 + 8. Rib centre surfaces and segments
    # Always compute faces if either visualisation is needed
    faces, edges = [], []
    if params.vis_rib_centre_surfaces or params.vis_rib_surface_segments:
        from slab_utils import create_rib_centre_surfaces
        faces, edges = create_rib_centre_surfaces(
            wing.Shape, all_center_lines, params.plane_normal
        )

    if params.vis_rib_centre_surfaces:
        from viz_utils import show_rib_centre_surfaces, show_rib_centre_edges
        if faces:
            show_rib_centre_surfaces(faces, doc, color=(0.8, 0.4, 0.8), transparency=50)
            print(f"Created {len(faces)} rib centre surfaces.")
        if edges:
            show_rib_centre_edges(edges, doc, line_color=(0.2, 0.5, 1.0), line_width=2)
            print(f"Created {len(edges)} rib centre edges.")
        if not faces and not edges:
            print("No rib centre geometry found.")

    if params.vis_rib_surface_segments:
        from slab_utils import split_rib_faces_by_crossings
        from viz_utils import show_rib_segments
        if faces:
            segments = split_rib_faces_by_crossings(faces)
            show_rib_segments(segments, doc)
            print(f"Created {len(segments)} rib face segments.")
        else:
            print("No faces to segment.")

    # holed_faces = add_ellipse_holes_to_faces(segments, margin=2.0)
    # show_rib_centre_surfaces(holed_faces, doc, color=(0.2,0.8,0.4), transparency=50)

    if True:
        from slab_utils import create_rectangular_cutout_from_boundary
        cutouts = []
        for seg in segments:
            cut = create_rectangular_cutout_from_boundary(seg)
            
            if cut:
                cutouts.append(cut)
                print("cut created")
            else:
                print("no cut")
        if cutouts:
            from viz_utils import show_rib_centre_surfaces
            show_rib_centre_surfaces(cutouts, doc, color=(1.0,0.5,0.0), transparency=30)
            print(f"Created {len(cutouts)} rectangular cutout faces.")

    fit_view(doc)
    doc.save()
    print("Document saved.")

