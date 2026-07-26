import FreeCAD
import FreeCADGui
import Part
import Mesh
from slab_utils import create_cut_result, merge_and_show_final, create_rib_centre_surfaces, split_rib_faces_by_crossings, cutouts_from_segmens, log_op
from point_utils import collect_rib_midpoints, create_rib_wires
from prism_utils import create_bridges_trimmed_to_wing
from viz_utils import fit_view
import MeshPart
import math


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
                 input_step_path="",
                 output_stl_path="",
                 doc_path=r"C:\Users\natha\git\ContinuousPath\wing.FCStd",
                 obj_name='Pad',
                 rib_spacing=5.0, xy_rib_width=0.1, rib_width=None, rib_angle=30.0,
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
                 vis_rect_cutouts=True):
        if construction_plane not in PLANE_DEFS:
            raise ValueError(f"construction_plane must be one of {list(PLANE_DEFS.keys())}")
        self.nozzle_diameter = nozzle_diameter
        self.wall_thickness = wall_thickness or nozzle_diameter
        self.rib_spacing = rib_spacing
        self.rib_angle = rib_angle
        self.doc_path = doc_path
        self.obj_name = obj_name
        self.input_step_path = input_step_path
        self.output_stl_path = output_stl_path
        self.grid_orientation = grid_orientation
        self.construction_plane = construction_plane
        pdef = PLANE_DEFS[construction_plane]
        self.plane_normal = pdef['normal']
        self.plane_axis_u = pdef['axis_u']
        self.plane_axis_v = pdef['axis_v']
        self.primary_dir = self._project_primary(primary_dir)
        self.z_step = z_step
        self.cutout_marging = cutout_marging
        self.vis_cut_wing = vis_cut_wing
        self.vis_centre_lines = vis_centre_lines
        self.vis_midpoints = vis_midpoints
        self.vis_wires = vis_wires
        self.vis_bridges = vis_bridges
        self.vis_rib_surface_segments = vis_rib_surface_segments
        self.vis_rib_centre_surfaces = vis_rib_centre_surfaces
        self.vis_final_solid = vis_final_solid
        self.vis_rect_cutouts = vis_rect_cutouts
        if xy_rib_width:
            self.rib_width = xy_rib_width / math.sin(math.radians(90-rib_angle))
        else:
            self.rib_width = rib_width

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
    # Load the wing geometry (STEP or FCStd)
    if params.input_step_path and (params.input_step_path.lower().endswith('.step') or params.input_step_path.lower().endswith('.stp')):
        doc = FreeCAD.newDocument("Wing")
        shape = Part.Shape()
        shape.read(params.input_step_path)
        wing_feature = doc.addObject("Part::Feature", "ImportedShape")
        wing_feature.Shape = shape
        doc.recompute()
        wing_shape = shape
    else:
        doc = FreeCAD.open(params.doc_path)
        wing = doc.getObject(params.obj_name)
        if not wing:
            raise RuntimeError("Object not found.")
        wing_shape = wing.Shape

    bb = wing_shape.BoundBox

    # 1. Compute grid lines and centre lines
    from slab_utils import create_angled_grid_lines
    log_op("Starting pipeline: Grid line generation")
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_center_lines = lines1 + lines2
    log_op(f"Generated {len(all_center_lines)} center lines")

    # 2. Generate rib centre surfaces and segments
    log_op("Generating rib centre surfaces")
    rib_centre_surfaces, _ = create_rib_centre_surfaces(wing_shape, all_center_lines, params.plane_normal, doc, vis=params.vis_rib_centre_surfaces)
    log_op("Splitting rib faces by crossings")
    rib_centre_segments = split_rib_faces_by_crossings(rib_centre_surfaces, doc, vis=params.vis_rib_surface_segments)

    # 3. Create cutout faces (rectangular holes)
    log_op("Creating cutout faces")
    cutout_faces = cutouts_from_segmens(rib_centre_segments, doc, vis=params.vis_rect_cutouts, margin=params.cutout_marging)

    # 4. Create cut result with holes (subtract holes from ribs)
    log_op("Creating cut result with holes")
    cut_result, all_center_lines, valid_ribs = create_cut_result(
        wing_shape, params, doc,
        vis_cut_wing=params.vis_cut_wing,
        vis_centre_lines=params.vis_centre_lines,
        cutout_faces=cutout_faces
    )

    # 5. Bridges (midpoints, wires, bridges)
    log_op("Collecting rib midpoints")
    points_by_rib = collect_rib_midpoints(wing_shape=wing_shape,
                                          rib_center_lines=all_center_lines,
                                          plane_normal=params.plane_normal,
                                          z_min=bb.ZMin + 0.5,
                                          z_max=bb.ZMax - 0.5,
                                          z_step=params.z_step,
                                          doc=doc,
                                          vis=params.vis_midpoints)
    log_op("Creating rib wires")
    wires = create_rib_wires(points_by_rib, doc, vis=params.vis_wires)
    log_op("Generating bridges")
    bridges = create_bridges_trimmed_to_wing(data_by_rib=points_by_rib,
                                             rib_center_lines=all_center_lines,
                                             plane_normal=params.plane_normal,
                                             rib_width=params.rib_width,
                                             bridge_height=0.5,
                                             wing_shape=wing_shape,
                                             extend_length=5.0,
                                             doc=doc,
                                             vis=params.vis_bridges)

    # 6. Final solid
    log_op("Merging final solid")
    final_solid = merge_and_show_final(cut_result, bridges, doc, vis=params.vis_final_solid)
    if final_solid is None or final_solid.isNull():
        print("Error: final solid is null – cannot proceed.")
        return
    print("Final solid (wing + bridges) created successfully.")

    # 7. Export to STL if output path is provided
    if params.output_stl_path:
        try:
            shape = final_solid
            # If the shape is a compound, fuse its components into a single solid
            if shape.ShapeType in ("Compound", "CompSolid"):
                solids = shape.Solids
                if solids:
                    fused = solids[0]
                    for s in solids[1:]:
                        fused = fused.fuse(s)
                    shape = fused
                else:
                    shape = shape.removeSplitter()
            # Convert shape to mesh using MeshPart (preferred) or fallback to temporary object
            try:
                mesh = MeshPart.meshFromShape(shape, LinearDeflection=0.1, AngularDeflection=0.5236)
                mesh.write(params.output_stl_path)
            except ImportError:
                # Fallback: create a temporary mesh object
                mesh_obj = FreeCAD.ActiveDocument.addObject("Mesh::Feature", "__temp_mesh__")
                mesh_obj.Mesh = Mesh.Mesh(shape)
                Mesh.export([mesh_obj], params.output_stl_path)
                FreeCAD.ActiveDocument.removeObject(mesh_obj.Name)
            print(f"Exported STL to {params.output_stl_path}")
        except Exception as e:
            print(f"Export failed: {e}")

    # 8. Save the FreeCAD document (optional)
    if params.input_step_path:
        if params.output_stl_path:
            fcstd_path = params.output_stl_path.replace('.stl', '.FCStd')
            doc.saveAs(fcstd_path)
            print(f"Saved FCStd document to {fcstd_path}")
    else:
        doc.save()
        print("Document saved.")

    fit_view(doc)


if __name__ == "__main__":
    params = LWInfillParams(
        # ===== INPUT / OUTPUT =====
        # input_step_path=r"C:\Users\natha\git\ContinuousPath\.in\test_wing.step",   # set to empty to use FCStd
        # output_stl_path=r"C:\Users\natha\Prints\test_wing\wing-FullGrid.stl",
        # Fallback for FCStd
        doc_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin.FCStd",
        # obj_name='Nose001',
        # obj_name='WingR1_msv001_solid',
        obj_name='Part__Feature_solid',
        # obj_name='Extrude001',
        # ===== GEOMETRY PARAMETERS =====
        nozzle_diameter=0.4,
        rib_spacing=20.0,
        xy_rib_width=0.13,
        rib_angle=30.0,
        grid_orientation=0.0,
        z_step=10,
        primary_dir=FreeCAD.Vector(0, 0, 1),
        construction_plane='XZ',
        # ===== VISUALISATION TOGGLES =====
        vis_cut_wing=False,
        vis_centre_lines=False,
        vis_midpoints=False,
        vis_wires=False,
        vis_bridges=False,
        vis_rib_centre_surfaces=False,
        vis_rib_surface_segments=False,
        vis_final_solid=True,
        vis_rect_cutouts=True,
        cutout_marging=0
    )
    main(params)