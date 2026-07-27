# main.py – fully optimised using trimesh for everything

import FreeCAD
import FreeCADGui
import Part
import Mesh
import math
import tempfile
import os
import numpy as np
import trimesh

# ---- Import only needed functions from our modules ----
from slab_utils import (
    create_angled_grid_lines,
    shape_to_trimesh,
    trimesh_to_freecad,
    get_rib_boundary_z_range,
    collect_midpoints,          # fast midpoint collection
    create_bridge_meshes,       # fast bridge creation
    create_holed_rib_union,     # for rib generation with holes
    repair_mesh,
)
from viz_utils import fit_view, show_mesh, show_rib_centre_lines

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
                 create_holes=True):
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
    # Load wing geometry (solid) and convert to trimesh
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

    # ---- Convert wing to trimesh once ----
    print("Converting wing to trimesh...")
    wing_mesh = shape_to_trimesh(wing_shape, deflection=0.5, angular=0.5)
    if wing_mesh is None:
        raise RuntimeError("Failed to convert wing to trimesh")
    print(f"Wing mesh: {len(wing_mesh.vertices)} vertices, {len(wing_mesh.faces)} faces")
    print(f"Is watertight: {wing_mesh.is_watertight}, is volume: {wing_mesh.is_volume}")

    # ---- 1. Generate grid lines (numpy) ----
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_lines = lines1 + lines2  # list of (start, end) numpy arrays
    print(f"Generated {len(all_lines)} center lines")

    # Convert to FreeCAD lines for visualisation (if needed)
    all_center_lines = []
    for start, end in all_lines:
        line = Part.makeLine(FreeCAD.Vector(*start), FreeCAD.Vector(*end))
        all_center_lines.append(line)

    # ---- 2. Create rib union (with optional holes) ----
    y_extent = (bb.YMax - bb.YMin) * 1.2
    rib_union = create_holed_rib_union(
        wing_mesh, all_lines, params,
        np.array(params.plane_normal),
        params.rib_width, y_extent,
        params.cutout_marging
    )
    if rib_union is None:
        print("No ribs generated – exiting.")
        return

    # ---- 3. Subtract ribs from wing ----
    print("Subtracting ribs from wing...")
    if wing_mesh.is_watertight and rib_union.is_volume:
        try:
            cut_wing = trimesh.boolean.difference([wing_mesh, rib_union], engine='manifold')
        except:
            try:
                cut_wing = trimesh.boolean.difference([wing_mesh, rib_union], engine='scad')
            except:
                cut_wing = wing_mesh
    else:
        cut_wing = wing_mesh

    # Visualise cut wing (if requested)
    cut_wing_fc = trimesh_to_freecad(cut_wing)
    if params.vis_cut_wing:
        show_mesh(cut_wing_fc, doc, "CutWingMesh", transparency=80)

    # ---- 4. Fast midpoints and bridges ----
    print("Computing midpoints for bridges...")
    # For each rib, get its vertical range inside the wing
    rib_z_ranges = []
    for start, end in all_lines:
        z_range = get_rib_boundary_z_range(wing_mesh, start, end, np.array(params.plane_normal))
        if z_range is None:
            rib_z_ranges.append((None, None))
        else:
            rib_z_ranges.append(z_range)

    # Collect midpoints along each rib using fast sampling
    data_by_rib = collect_midpoints(all_lines, rib_z_ranges, params.z_step)

    # Build bridge meshes
    bridge_meshes = create_bridge_meshes(
        data_by_rib,
        np.array(params.plane_normal),
        params.rib_width,
        bridge_height=0.5
    )
    if bridge_meshes:
        print(f"Unioning {len(bridge_meshes)} bridges...")
        bridge_union = bridge_meshes[0]
        for bm in bridge_meshes[1:]:
            try:
                bridge_union = trimesh.boolean.union([bridge_union, bm], engine='manifold')
            except:
                try:
                    bridge_union = trimesh.boolean.union([bridge_union, bm], engine='scad')
                except:
                    print("  Skipping a bridge union due to failure")
                    continue
        # Trim bridges to wing (optional, but recommended)
        print("Trimming bridges to wing...")
        try:
            bridge_union = bridge_union.intersection(wing_mesh, engine='manifold')
        except:
            try:
                bridge_union = bridge_union.intersection(wing_mesh, engine='scad')
            except:
                pass  # keep full bridge
    else:
        bridge_union = None

    # ---- 5. Final merge ----
    if bridge_union is not None:
        print("Merging wing with bridges...")
        try:
            final_mesh = trimesh.boolean.union([cut_wing, bridge_union], engine='manifold')
        except:
            try:
                final_mesh = trimesh.boolean.union([cut_wing, bridge_union], engine='scad')
            except:
                final_mesh = cut_wing
    else:
        final_mesh = cut_wing

    # ---- 6. Visualise final ----
    final_fc = trimesh_to_freecad(final_mesh)
    if params.vis_final_solid:
        show_mesh(final_fc, doc, "FinalMesh", color=(0.8,0.8,0.8), transparency=50)

    # ---- 7. Export STL ----
    if params.output_stl_path:
        final_mesh.export(params.output_stl_path)
        print(f"Exported STL to {params.output_stl_path}")

    # ---- 8. Save FCStd ----
    if params.input_step_path:
        if params.output_stl_path:
            fcstd_path = params.output_stl_path.replace('.stl', '.FCStd')
            doc.saveAs(fcstd_path)
            print(f"Saved FCStd document to {fcstd_path}")
    else:
        doc.save()
        print("Document saved.")

    fit_view(doc)

# ------------------------------------------------------------
# Run with your parameters
# ------------------------------------------------------------
if __name__ == "__main__":
    params = LWInfillParams(
        doc_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin.FCStd",
        # obj_name='Part__Feature_solid',
        obj_name='WingR1_msv_orient001_solid',
        nozzle_diameter=0.4,
        rib_spacing=20.0,
        xy_rib_width=0.13,
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
        cutout_marging=0,
        create_holes=False   # set to True if you want cutouts
    )
    main(params)