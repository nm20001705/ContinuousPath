# main.py – fully optimised using trimesh for everything

import FreeCAD
import Part
import math
import numpy as np

# ---- Import only needed functions from our modules ----
# ---- Import the new function ----
from slab_utils import (
    create_rib_surfaces_trimesh,
    shape_to_trimesh, 
    create_angled_grid_lines,
    clip_surfaces_to_solid, 
    build_rib_segments_analytical
)

from bridge_utils import create_bridges_analytical

from viz_utils import fit_view, show_rib_centre_lines

# ---- PLANE_DEFS and LWInfillParams (with create_holes) ----
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
                 input_step_path="", output_stl_path="",
                 doc_path=r"C:\Users\natha\git\ContinuousPath\wing.FCStd",
                 obj_name='Pad',
                 rib_spacing=5.0, xy_rib_width=0.1, rib_width=None, rib_angle=30.0,
                 grid_orientation=0.0, primary_dir=None, z_step=10, cutout_marging=2,
                 construction_plane='XZ',
                 vis_cut_wing=True, vis_rib_centre_surfaces=True,
                 vis_centre_lines=True, vis_midpoints=True,
                 vis_wires=True, vis_bridges=True,
                 vis_rib_surface_segments=True,
                 vis_final_solid=True, vis_rect_cutouts=True,
                 create_holes=True, 
                 vis_rib_centre_surfaces_clip=True, 
                 vis_rib_segments=True, 
                 bridge_height=0.4, 
                 vis_bridge=True):
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
        self.create_holes = create_holes
        self.vis_rib_centre_surfaces_clip=vis_rib_centre_surfaces_clip
        self.vis_rib_segments=vis_rib_segments
        self.bridge_height=bridge_height
        self.vis_bridge=vis_bridge
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

# ------------------------------------------------------------
# Main function – fully optimised
# ------------------------------------------------------------
def main(params):
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
        wing_obj = doc.getObject(params.obj_name)
        if not wing_obj:
            raise RuntimeError("Object not found.")
        wing_shape = wing_obj.Shape

    bb = wing_shape.BoundBox

    # ---- Convert wing to trimesh (for later use, e.g., boolean ops) ----
    print("Converting wing to trimesh...")
    wing_mesh = shape_to_trimesh(wing_shape, deflection=0.5, angular=0.5)
    if wing_mesh is None:
        raise RuntimeError("Failed to convert wing to trimesh")
    print(f"Wing bounds: {wing_mesh.bounds}")
    print(f"Is watertight: {wing_mesh.is_watertight}, is volume: {wing_mesh.is_volume}")

    # Test section at center of bounding box
    center = wing_mesh.bounds.mean(axis=0)
    test_section = wing_mesh.section(plane_origin=center, plane_normal=[0,0,1])
    print(f"Test section at center {center}: {test_section is not None}")
    # ---- Generate grid lines (numpy arrays) ----
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_lines_np = lines1 + lines2   # list of (start, end) numpy arrays
    print(f"Generated {len(all_lines_np)} rib center lines")

    # Convert to FreeCAD lines for the sectioning function
    all_center_lines_fc = []
    for start, end in all_lines_np:
        line = Part.makeLine(FreeCAD.Vector(*start), FreeCAD.Vector(*end))
        all_center_lines_fc.append(line)

    # ---- Create rib centre surfaces using fast FreeCAD section + meshing ----
    plane_normal_np = np.array([params.plane_normal.x, params.plane_normal.y, params.plane_normal.z])
    rib_faces = create_rib_surfaces_trimesh(
        wing_mesh,
        all_lines_np,
        plane_normal_np,
        doc=doc,
        vis=params.vis_rib_centre_surfaces
    )
    print(f"Created {len(rib_faces)} rib centre faces")

    # ---- Get the exact rib cross‑sections ----
    clipped_ribs, rib_indices = clip_surfaces_to_solid(
        all_lines_np,               # rib centerlines (numpy arrays)
        wing_mesh,                  # the wing mesh
        plane_normal_np,            # construction plane normal
        doc=doc,
        vis=params.vis_rib_centre_surfaces_clip
    )

    active_lines = [all_lines_np[i] for i in rib_indices]

    # Now split
    rib_segments = build_rib_segments_analytical(
        wing_mesh,
        all_lines_np,
        plane_normal_np,
        doc=doc,
        vis=params.vis_rib_segments
    )
    primary_dir_np = np.array([
        params.primary_dir.x,
        params.primary_dir.y,
        params.primary_dir.z
    ])
    primary_dir_np /= np.linalg.norm(primary_dir_np)
    
    # Build bridges
    bridge_mesh = create_bridges_analytical(
        wing_mesh,
        all_lines_np,
        primary_dir_np,
        z_step=params.z_step,                       # your chosen step
        bridge_height=params.bridge_height,
        construction_plane_normal=plane_normal_np,   # <-- ADD THIS
        doc=doc,
        vis=True
    )

    # ---- Show centre lines if requested ----
    if params.vis_centre_lines:
        show_rib_centre_lines(all_center_lines_fc, doc)


    # ---- Save document ----
    doc.save()
    print("Document saved.")
    fit_view(doc)
# ------------------------------------------------------------
# Run with your parameters
# ------------------------------------------------------------
if __name__ == "__main__":
    params = LWInfillParams(
        doc_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin.FCStd",
        obj_name='Part__Feature_solid',
        # obj_name='WingR1_msv_orient001_solid',
        nozzle_diameter=0.4,
        rib_spacing=20.0,
        bridge_height = 0.4,
        xy_rib_width=0.13,
        rib_angle=30.0,
        grid_orientation=0.0,
        z_step=1,
        primary_dir=FreeCAD.Vector(0, 0, 1),
        construction_plane='XZ',
        vis_cut_wing=False,
        vis_centre_lines=False,
        vis_midpoints=False,
        vis_wires=False,
        vis_bridges=False,
        vis_rib_centre_surfaces=False,
        vis_rib_centre_surfaces_clip=False,
        vis_rib_surface_segments=False,
        vis_final_solid=False,
        vis_rect_cutouts=True,
        cutout_marging=0,
        create_holes=False,    # set to True if you want cutouts
        vis_rib_segments = True, 
        vis_bridge = True
    )
    main(params)
