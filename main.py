# main.py – cleaned version with only used parameters

from types import SimpleNamespace
import FreeCAD
import Part
import numpy as np

from slab_utils import (
    create_rib_surfaces_trimesh,
    shape_to_trimesh,
    create_angled_grid_lines,
    clip_surfaces_to_solid,
    build_rib_segments_analytical,
)
from bridge_utils import create_bridges_analytical
from viz_utils import fit_view, show_rib_centre_lines
from hole_utils import create_holes

# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------
def main(params):
    # Load or open the document
    if params.input_step_path and params.input_step_path.lower().endswith(('.step', '.stp')):
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

    # Convert wing to trimesh
    print("Converting wing to trimesh...")
    wing_mesh = shape_to_trimesh(wing_shape, deflection=0.5, angular=0.5)
    if wing_mesh is None:
        raise RuntimeError("Failed to convert wing to trimesh")
    print(f"Wing bounds: {wing_mesh.bounds}")
    print(f"Is watertight: {wing_mesh.is_watertight}, is volume: {wing_mesh.is_volume}")

    # Test section (can be removed later)
    center = wing_mesh.bounds.mean(axis=0)
    test_section = wing_mesh.section(plane_origin=center, plane_normal=[0, 0, 1])
    print(f"Test section at center {center}: {test_section is not None}")

    # Generate rib centre lines
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_lines_np = lines1 + lines2
    print(f"Generated {len(all_lines_np)} rib center lines")

    # Convert to FreeCAD lines (for optional visualization)
    all_center_lines_fc = []
    for start, end in all_lines_np:
        line = Part.makeLine(FreeCAD.Vector(*start), FreeCAD.Vector(*end))
        all_center_lines_fc.append(line)

    # Numpy versions of plane vectors
    plane_normal_np = np.array([params.plane_normal.x, params.plane_normal.y, params.plane_normal.z])
    primary_dir_np = np.array([params.primary_dir.x, params.primary_dir.y, params.primary_dir.z])
    primary_dir_np /= np.linalg.norm(primary_dir_np)

    # Rib centre surfaces (visualization only)
    rib_faces = create_rib_surfaces_trimesh(
        wing_mesh,
        all_lines_np,
        plane_normal_np,
        doc=doc,
        vis=params.vis_rib_centre_surfaces,
    )
    print(f"Created {len(rib_faces)} rib centre faces")

    # Clipped rib surfaces
    clipped_ribs, rib_indices = clip_surfaces_to_solid(
        all_lines_np,
        wing_mesh,
        plane_normal_np,
        doc=doc,
        vis=params.vis_rib_centre_surfaces_clip,
    )

    # Rib segments (the actual pieces inside the wing)
    rib_segments, segment_bounds = build_rib_segments_analytical(
        wing_mesh,
        all_lines_np,
        plane_normal_np,
        doc=doc,
        vis=params.vis_rib_segments,
    )

    # Precompute slab normals for every rib line
    lines_data = []
    for start, end in all_lines_np:
        d = end - start
        if np.linalg.norm(d) < 1e-8:
            lines_data.append(None)
            continue
        dir_rib = d / np.linalg.norm(d)
        slab_normal = np.cross(dir_rib, plane_normal_np)
        slab_normal /= np.linalg.norm(slab_normal)
        lines_data.append({'point': start, 'dir': dir_rib, 'slab_normal': slab_normal})

    # Convert segment bounds into the format needed by the bridge/hole functions
    bridge_segments = []
    for line_idx, s0, s1 in segment_bounds:
        ld = lines_data[line_idx]
        p0 = ld['point'] + s0 * ld['dir']
        p1 = ld['point'] + s1 * ld['dir']
        bridge_segments.append({
            'p0': p0,
            'p1': p1,
            'dir': ld['dir'],
            'slab_normal': ld['slab_normal'],
        })

    # Create bridges
    bridge_mesh = create_bridges_analytical(
        wing_mesh,
        bridge_segments,
        primary_dir_np,
        z_step=1,
        bridge_height=params.bridge_height,
        doc=doc,
        vis=params.vis_bridge,
    )

    # Create holes
    def hole_condition(x):
        return np.sqrt(max(0, 1 - (2 * x - 1) ** 2))

    hole_mesh = create_holes(
        wing_mesh,
        bridge_segments,
        primary_dir_np,
        z_step=1,
        point_condition=hole_condition,
        doc=doc,
        vis=params.vis_hole,
        hole_margin=params.hole_margin,
    )

    # Show centre lines if desired
    if params.vis_centre_lines:
        show_rib_centre_lines(all_center_lines_fc, doc)

    doc.save()
    print("Document saved.")
    fit_view(doc)


# ------------------------------------------------------------
# Runtime configuration – change your values here
# ------------------------------------------------------------
if __name__ == "__main__":
    # Construction plane definitions (needed to derive normal & axes)
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

    # Choose your plane
    construction_plane = 'XZ'
    pdef = PLANE_DEFS[construction_plane]

    # Base parameters
    params = SimpleNamespace(
        doc_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin.FCStd",
        obj_name='Part__Feature_solid',
        rib_spacing=20.0,
        rib_angle=30.0,
        grid_orientation=0.0,
        primary_dir=FreeCAD.Vector(0, 0, 1),
        bridge_height=0.4,
        hole_margin=0.5,
        input_step_path="",          # set to a .step/.stp file if needed
        # Visibility toggles
        vis_rib_centre_surfaces=False,
        vis_rib_centre_surfaces_clip=False,
        vis_rib_segments=True,
        vis_centre_lines=False,
        vis_bridge=True,
        vis_hole=True,
    )

    # Attach plane vectors to the namespace
    params.plane_normal = pdef['normal']
    params.plane_axis_u = pdef['axis_u']
    params.plane_axis_v = pdef['axis_v']

    # Project primary_dir onto the construction plane if needed
    if params.primary_dir is not None:
        n = params.plane_normal
        pd = params.primary_dir
        dot = pd.x * n.x + pd.y * n.y + pd.z * n.z
        proj = FreeCAD.Vector(pd.x - dot * n.x, pd.y - dot * n.y, pd.z - dot * n.z)
        if proj.Length > 1e-6:
            params.primary_dir = proj.normalize()
        else:
            params.primary_dir = None

    main(params)