# main.py – fully optimised with cached wing slices

from types import SimpleNamespace
import FreeCAD
import Part
import math
import numpy as np
import trimesh
from shapely.geometry import LineString, Polygon as ShapelyPolygon
from shapely.ops import linemerge, polygonize
from scipy.spatial import ConvexHull

from slab_utils import (
    create_rib_surfaces_trimesh,
    shape_to_trimesh,
    create_angled_grid_lines,
    clip_surfaces_to_solid,
    build_rib_segments_analytical,
)
from bridge_utils import create_bridges_analytical
from hole_utils import create_holes
from viz_utils import fit_view, show_rib_centre_lines


# ------------------------------------------------------------
# Precompute wing cross‑sections for all Z‑slices
# ------------------------------------------------------------
def precompute_slices(wing_mesh, prim, z_step, d_min, d_max):
    """
    Returns
    -------
    z_vals : np.ndarray   – sorted Z coordinates of the slices
    slices : list         – list of (polygon_2d, to_3d_matrix) for each Z
    """
    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)

    z_vals = []
    slices = []

    d = d_min
    while d <= d_max + 1e-9:
        plane_origin = d * prim
        try:
            result = trimesh.intersections.mesh_plane(
                wing_mesh, plane_normal=prim, plane_origin=plane_origin
            )
        except Exception:
            d += z_step
            continue
        if isinstance(result, tuple):
            lines = result[0]
        else:
            lines = result
        if lines is None or len(lines) == 0:
            d += z_step
            continue

        origin = np.asarray(plane_origin)
        segments_2d = []
        for seg in lines:
            p1 = seg[0] - origin
            p2 = seg[1] - origin
            u1, v1 = np.dot(p1, u_ax), np.dot(p1, v_ax)
            u2, v2 = np.dot(p2, u_ax), np.dot(p2, v_ax)
            if np.linalg.norm([u2 - u1, v2 - v1]) > 1e-9:
                segments_2d.append(LineString([(u1, v1), (u2, v2)]))

        if not segments_2d:
            d += z_step
            continue
        merged = linemerge(segments_2d)
        if merged.is_empty:
            d += z_step
            continue
        polys = list(polygonize(merged))
        if not polys:
            all_pts = []
            for line in segments_2d:
                all_pts.extend(line.coords)
            pts = np.array(all_pts)
            if len(pts) < 3:
                d += z_step
                continue
            try:
                hull = ConvexHull(pts)
                poly = ShapelyPolygon(pts[hull.vertices])
            except Exception:
                d += z_step
                continue
        else:
            poly = max(polys, key=lambda p: p.area)

        if poly.is_empty or poly.area < 1e-8:
            d += z_step
            continue

        to_3d = np.eye(4)
        to_3d[:3, 0] = u_ax
        to_3d[:3, 1] = v_ax
        to_3d[:3, 3] = origin

        z_vals.append(d)
        slices.append((poly, to_3d))
        d += z_step

    return np.array(z_vals), slices

# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------
def main(params):
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

    # ---- Convert wing to trimesh ----
    print("Converting wing to trimesh...")
    wing_mesh = shape_to_trimesh(wing_shape, deflection=0.5, angular=0.5)
    if wing_mesh is None:
        raise RuntimeError("Failed to convert wing to trimesh")
    print(f"Wing bounds: {wing_mesh.bounds}")
    print(f"Is watertight: {wing_mesh.is_watertight}, is volume: {wing_mesh.is_volume}")

    # ---- Plane vectors (used everywhere) ----
    plane_normal_np = np.array([params.plane_normal.x, params.plane_normal.y, params.plane_normal.z])
    primary_dir_np = np.array([params.primary_dir.x, params.primary_dir.y, params.primary_dir.z])
    primary_dir_np /= np.linalg.norm(primary_dir_np)

    # ---- Slice range (Z) ----
    bbox = wing_mesh.bounds
    min_pt, max_pt = bbox[0], bbox[1]
    corners = np.array([
        [min_pt[0], min_pt[1], min_pt[2]],
        [min_pt[0], min_pt[1], max_pt[2]],
        [min_pt[0], max_pt[1], min_pt[2]],
        [min_pt[0], max_pt[1], max_pt[2]],
        [max_pt[0], min_pt[1], min_pt[2]],
        [max_pt[0], min_pt[1], max_pt[2]],
        [max_pt[0], max_pt[1], min_pt[2]],
        [max_pt[0], max_pt[1], max_pt[2]],
    ])
    proj = np.dot(corners, primary_dir_np)
    d_min = proj.min() - 1e-3
    d_max = proj.max() + 1e-3

    # ---- Precompute all Z-slices (once!) ----
    print("Precomputing wing slices...")
    z_vals, slices = precompute_slices(wing_mesh, primary_dir_np, params.z_step, d_min=d_min, d_max=d_max)
    print(f"Precomputed {len(slices)} slices")

    # ---- Generate rib centre lines ----
    bb = wing_shape.BoundBox
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_lines_np = lines1 + lines2
    print(f"Generated {len(all_lines_np)} rib center lines")

    # Convert to FreeCAD lines (optional visualization)
    all_center_lines_fc = []
    for start, end in all_lines_np:
        line = Part.makeLine(FreeCAD.Vector(*start), FreeCAD.Vector(*end))
        all_center_lines_fc.append(line)

    # ---- Rib centre surfaces (visualization only) ----
    rib_faces = create_rib_surfaces_trimesh(
        wing_mesh,
        all_lines_np,
        plane_normal_np,
        doc=doc,
        vis=params.vis_rib_centre_surfaces,
    )
    print(f"Created {len(rib_faces)} rib centre faces")

    # ---- Clipped rib surfaces (if needed for vis) ----
    if params.vis_rib_centre_surfaces_clip:
        clipped_ribs, rib_indices = clip_surfaces_to_solid(
            all_lines_np,
            wing_mesh,
            plane_normal_np,
            doc=doc,
            vis=True,
        )

    # ---- Build rib segments (the actual pieces inside the wing) ----
    rib_segments, segment_bounds = build_rib_segments_analytical(
        wing_mesh,
        all_lines_np,
        plane_normal_np,
        doc=doc,
        vis=params.vis_rib_segments,
    )

    # ---- Precompute slab normals for all lines ----
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

    # ---- Convert segment bounds into bridge/hole format ----
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

    # ---- Create bridges ----
    bridge_mesh = create_bridges_analytical(
        wing_mesh,
        bridge_segments,
        primary_dir_np,
        z_step=1,
        bridge_height=params.bridge_height,
        z_vals=z_vals,          # <-- new parameter
        slices=slices,          # list, not dict
        doc=doc,
        vis=params.vis_bridge,
    )

    # ---- Create holes ----
    def hole_condition(x):
        return np.sqrt(max(0, 1 - (2 * x - 1) ** 2))

    hole_mesh = create_holes(
        wing_mesh,
        bridge_segments,
        primary_dir_np,
        z_step=1,
        point_condition=hole_condition,
        z_vals=z_vals,          # <-- new parameter
        slices=slices,
        doc=doc,
        vis=params.vis_hole,
        hole_margin=params.hole_margin,
    )

    # ---- Show centre lines if requested ----
    if params.vis_centre_lines:
        show_rib_centre_lines(all_center_lines_fc, doc)

    # ---- Save ----
    doc.save()
    print("Document saved.")
    fit_view(doc)


# ------------------------------------------------------------
# Runtime configuration
# ------------------------------------------------------------
if __name__ == "__main__":
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

    construction_plane = 'XZ'
    pdef = PLANE_DEFS[construction_plane]

    params = SimpleNamespace(
        doc_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin.FCStd",
        obj_name='Part__Feature_solid',
        rib_spacing=20.0,
        rib_angle=30.0,
        grid_orientation=0.0,
        primary_dir=FreeCAD.Vector(0, 0, 1),
        bridge_height=0.4,
        hole_margin=0.5,
        input_step_path="",
        vis_rib_centre_surfaces=False,
        vis_rib_centre_surfaces_clip=False,
        vis_rib_segments=True,
        vis_centre_lines=False,
        vis_bridge=True,
        vis_hole=True,
        z_step=0.2
    )

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