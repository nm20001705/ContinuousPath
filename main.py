import FreeCAD
import FreeCADGui
import Part
from slab_utils import create_cut_result, merge_and_show_final, create_rectangular_cutout_from_boundary, create_rib_centre_surfaces, split_rib_faces_by_crossings, cutouts_from_segmens
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
    def __init__(self, nozzle_diameter=0.4, wall_thickness=None, doc_path=r"C:\Users\natha\git\ContinuousPath\wing.FCStd", obj_name='Pad',
                 rib_spacing=5.0, rib_width=None, rib_angle=30.0,
                 grid_orientation=0.0, primary_dir=None, z_step=10, cutout_marging=2,
                 construction_plane='XZ',
                 vis_cut_wing=True,
                 vis_rib_centre_surfaces=True,
                 vis_centre_lines=True,
                 vis_midpoints=True,
                 vis_wires=True,
                 vis_bridges=True,
                 vis_rib_surface_segments=True,
                 vis_final_solid=True, 
                 vis_rect_cutouts = True):
        if construction_plane not in PLANE_DEFS:
            raise ValueError(f"construction_plane must be one of {list(PLANE_DEFS.keys())}")
        self.nozzle_diameter = nozzle_diameter
        self.wall_thickness = wall_thickness or nozzle_diameter
        self.rib_spacing = rib_spacing
        self.rib_width = rib_width or (2.0 * nozzle_diameter)
        self.rib_angle = rib_angle
        self.doc_path = doc_path
        self.obj_name =obj_name
        self.grid_orientation = grid_orientation
        self.construction_plane = construction_plane
        pdef = PLANE_DEFS[construction_plane]
        self.plane_normal = pdef['normal']
        self.plane_axis_u = pdef['axis_u']
        self.plane_axis_v = pdef['axis_v']
        self.primary_dir = self._project_primary(primary_dir)
        self.z_step = z_step
        self.cutout_marging=cutout_marging
        self.vis_cut_wing = vis_cut_wing
        self.vis_centre_lines = vis_centre_lines
        self.vis_midpoints = vis_midpoints
        self.vis_wires = vis_wires
        self.vis_bridges = vis_bridges
        self.vis_rib_surface_segments = vis_rib_surface_segments
        self.vis_rib_centre_surfaces = vis_rib_centre_surfaces
        self.vis_final_solid = vis_final_solid
        self.vis_rect_cutouts = vis_rect_cutouts

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

def main(params):
    doc = FreeCAD.open(params.doc_path)
    wing = doc.getObject(params.obj_name)
    if not wing:
        raise RuntimeError("Object not found.")

    bb = wing.Shape.BoundBox

    # 1. Compute grid lines and centre lines
    from slab_utils import create_angled_grid_lines
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_center_lines = lines1 + lines2

    # 2. Generate rib centre surfaces and segments
    rib_centre_surfaces, _ = create_rib_centre_surfaces(wing.Shape, all_center_lines, params.plane_normal, doc, vis=params.vis_rib_centre_surfaces)
    rib_centre_segments = split_rib_faces_by_crossings(rib_centre_surfaces, doc, vis=params.vis_rib_surface_segments)

    # 3. Create cutout faces (rectangular holes)
    cutout_faces = cutouts_from_segmens(rib_centre_segments, doc, vis=params.vis_rect_cutouts, margin=params.cutout_marging)

    # 4. Create cut result with holes (the cutout faces are used to subtract holes from ribs)
    cut_result, all_center_lines, valid_ribs = create_cut_result(
        wing.Shape, params, doc,
        vis_cut_wing=params.vis_cut_wing,
        vis_centre_lines=params.vis_centre_lines,
        cutout_faces=cutout_faces
    )

    # 5. Bridges (midpoints, wires, bridges)
    points_by_rib = collect_rib_midpoints(wing_shape=wing.Shape, rib_center_lines=all_center_lines,
                                          plane_normal=params.plane_normal,
                                          z_min=bb.ZMin+0.5, z_max=bb.ZMax-0.5,
                                          z_step=params.z_step, doc=doc, vis=params.vis_midpoints)
    wires = create_rib_wires(points_by_rib, doc, vis=params.vis_wires)
    bridges = create_bridges_trimmed_to_wing(data_by_rib=points_by_rib, rib_center_lines=all_center_lines,
                                             plane_normal=params.plane_normal, rib_width=params.rib_width,
                                             bridge_height=0.5, wing_shape=wing.Shape, extend_length=5.0,
                                             doc=doc, vis=params.vis_bridges)

    # 6. Final solid
    final_solid = merge_and_show_final(cut_result, bridges, doc, vis=params.vis_final_solid)
    print("Final solid (wing + bridges) created successfully.")

    fit_view(doc)
    doc.save()
    print("Document saved.")

if __name__ == "__main__":
    params = LWInfillParams(
        obj_name='Scale',
        nozzle_diameter=0.4,
        rib_spacing=30.0,
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
        vis_rib_centre_surfaces=False,
        vis_rib_surface_segments=False,
        vis_final_solid=True,
        vis_rect_cutouts=True, 
        cutout_marging = 2
    )
    main(params=params)