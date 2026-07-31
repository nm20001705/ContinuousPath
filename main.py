# main.py – fully optimised with cached wing slices, mesh-hash protection, and robust guards

from types import SimpleNamespace
import FreeCAD
import Part
import math
import numpy as np
import trimesh
from shapely.geometry import LineString, Polygon as ShapelyPolygon
from shapely.ops import linemerge, polygonize
from scipy.spatial import ConvexHull
from tqdm import tqdm
import pickle
import os

from slab_utils import (
    create_rib_surfaces_trimesh,
    shape_to_trimesh,
    create_angled_grid_lines,
    clip_surfaces_to_solid,
    build_rib_segments_analytical,
)
from bridge_utils import create_bridges_analytical
from hole_utils import create_holes_analytical
from viz_utils import fit_view, show_rib_centre_lines


# ------------------------------------------------------------
# Precompute wing cross‑sections for all Z‑slices
# ------------------------------------------------------------
def precompute_slices(wing_mesh, prim, z_step, d_min, d_max):
    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)

    z_vals = []
    slices = []          # list of lists: each element is a list of (poly, to_3d)
    num_steps = int((d_max - d_min) / z_step) + 1
    pbar = tqdm(total=num_steps, desc="Precomputing slices")
    d = d_min
    while d <= d_max + 1e-9:
        plane_origin = d * prim
        try:
            result = trimesh.intersections.mesh_plane(
                wing_mesh, plane_normal=prim, plane_origin=plane_origin
            )
        except Exception:
            d += z_step
            pbar.update(1)
            continue
        if isinstance(result, tuple):
            lines = result[0]
        else:
            lines = result
        if lines is None or len(lines) == 0:
            d += z_step
            pbar.update(1)
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
            pbar.update(1)
            continue
        merged = linemerge(segments_2d)
        if merged.is_empty:
            d += z_step
            pbar.update(1)
            continue
        polys = list(polygonize(merged))
        if not polys:
            # Fallback: convex hull of all endpoints
            all_pts = []
            for line in segments_2d:
                all_pts.extend(line.coords)
            pts = np.array(all_pts)
            if len(pts) < 3:
                d += z_step
                pbar.update(1)
                continue
            try:
                hull = ConvexHull(pts)
                polys = [ShapelyPolygon(pts[hull.vertices])]
            except Exception:
                d += z_step
                pbar.update(1)
                continue

        # Build one to_3d matrix for this slice (shared by all polygons)
        to_3d = np.eye(4)
        to_3d[:3, 0] = u_ax
        to_3d[:3, 1] = v_ax
        to_3d[:3, 3] = origin

        # Store all non‑empty polygons with the shared matrix
        slice_polys = []
        for poly in polys:
            if poly.is_empty or poly.area < 1e-8:
                continue
            slice_polys.append(poly)

        if slice_polys:          # only keep slices that actually contain something
            z_vals.append(d)
            slices.append((slice_polys, to_3d))

        d += z_step
        pbar.update(1)

    pbar.close()
    return np.array(z_vals), slices

# ------------------------------------------------------------
# Improved mesh repair (preserves intentional holes)
# ------------------------------------------------------------
def repair_mesh(mesh, max_hole=3.0):
    """
    Minimal repair that fixes defects up to `max_hole` mm while preserving
    large intentional holes (servo cutouts, lightening holes).
    """
    if mesh is None:
        return None
    if isinstance(mesh, trimesh.Scene):
        meshes = [g for g in mesh.geometry.values() if isinstance(g, trimesh.Trimesh)]
        if not meshes:
            return None
        mesh = trimesh.util.concatenate(meshes) if len(meshes) > 1 else meshes[0]
        if mesh is None:
            return None
    if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0:
        return None

    # Clean up vertices and degenerate faces
    try:
        mesh.merge_vertices()
    except:
        pass
    try:
        mesh.remove_degenerate_faces()
    except:
        pass
    try:
        mesh.fix_normals()
    except:
        pass

    # Fill holes up to `max_hole` mm – close small defects but keep the servo hole
    try:
        mesh = trimesh.repair.fill_holes(mesh, max_hole=max_hole)
    except:
        pass
    try:
        mesh.fix_normals()
    except:
        pass

    return mesh

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

    # ---- Convert wing to trimesh (with minimal repair) ----
    print("Converting wing to trimesh...")
    wing_mesh = shape_to_trimesh(wing_shape, deflection=0.5, angular=0.5)
    if wing_mesh is None:
        raise RuntimeError("Failed to convert wing to trimesh")

    # ---- Mesh diagnostics ----
    print("=== Mesh diagnostics ===")
    print(f"Vertices: {len(wing_mesh.vertices)}")
    print(f"Faces: {len(wing_mesh.faces)}")
    print(f"Watertight: {wing_mesh.is_watertight}")
    print(f"Volume: {wing_mesh.is_volume}")

    # Number of disconnected pieces
    try:
        pieces = wing_mesh.split()
        if isinstance(pieces, list):
            print(f"Disconnected components: {len(pieces)}")
        else:
            print("Disconnected components: 1 (no split)")
    except Exception as e:
        print(f"Could not split mesh: {e}")

    print(f"Bounds: {wing_mesh.bounds}")

    # Non‑manifold edges (edges shared by ≠2 faces)
    try:
        edges_sorted = wing_mesh.edges_sorted
        unique_edges, counts = np.unique(edges_sorted, return_counts=True, axis=0)
        bad_edges = unique_edges[counts != 2]
        print(f"Non‑manifold edges: {len(bad_edges)}")
    except Exception as e:
        print(f"Could not count non‑manifold edges: {e}")

    # Face area statistics
    try:
        areas = wing_mesh.area_faces
        print(f"Face area range: {areas.min():.6f} to {areas.max():.3f}")
    except Exception as e:
        print(f"Could not compute face areas: {e}")

    # Boundary edges
    try:
        boundary_edges = wing_mesh.edges[wing_mesh.edges_unique, :]
        print(f"Boundary edges: {len(boundary_edges)}")
    except Exception as e:
        print(f"Could not compute boundary edges: {e}")

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

    # ---- Mesh hash for cache validation ----
    mesh_hash = hash(wing_mesh.vertices.tobytes())

    # ---- Precompute / load Z-slices from cache ----
    cache_dir = os.path.dirname(params.doc_path) if params.doc_path else os.getcwd()
    cache_filename = f"{params.obj_name}.pkl"
    cache_path = os.path.join(cache_dir, cache_filename)

    z_vals, slices = None, None
    if os.path.exists(cache_path):
        print(f"Loading cached slices from {cache_path}...")
        try:
            with open(cache_path, 'rb') as f:
                data = pickle.load(f)
            if (data.get('z_step') == params.z_step and
                data.get('mesh_hash') == mesh_hash):
                z_vals = data['z_vals']
                slices = data['slices']
                print(f"Loaded {len(slices)} slices from cache.")
            else:
                print("Cache is stale (mesh or z_step changed); recomputing...")
        except Exception as e:
            print(f"Failed to load cache ({e}); recomputing...")

    if z_vals is None:
        print("Precomputing wing slices...")
        z_vals, slices = precompute_slices(wing_mesh, primary_dir_np, params.z_step,
                                           d_min=d_min, d_max=d_max)
        print(f"Precomputed {len(slices)} slices")
        try:
            with open(cache_path, 'wb') as f:
                pickle.dump({'z_vals': z_vals, 'slices': slices,
                             'z_step': params.z_step, 'mesh_hash': mesh_hash}, f)
            print(f"Slices cached to {cache_path}")
        except Exception as e:
            print(f"Warning: could not save slice cache: {e}")

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
    bridge_mesh, bridge_solid = create_bridges_analytical(
        rib_segments,
        bridge_segments,
        primary_dir_np,
        z_vals=z_vals,
        bridge_height=params.bridge_height,
        margin=0.0,          # or params.bridge_margin if you want one
        doc=doc,
        vis=params.vis_bridge,
        thickness=params.thickness
    )

    # ---- Create holes ----
    def hole_condition(x):
        return np.sqrt(max(0, 1 - (2 * x - 1) ** 2))

    # def hole_condition(x):
    #     return 1

    hole_mesh, hole_solid = create_holes_analytical(
        rib_segments,          # <-- the trimesh list from build_rib_segments_analytical
        bridge_segments,       # <-- unchanged, same list you already build
        primary_dir_np,
        z_vals=z_vals,
        point_condition=hole_condition,
        hole_margin=params.hole_margin,
        doc=doc,
        vis=params.vis_hole,
        thickness=params.thickness
    )

    # ---- Show centre lines if requested ----
    if params.vis_centre_lines:
        show_rib_centre_lines(all_center_lines_fc, doc)

    # ---- Show wing mesh (for debugging) ----
    if params.vis_wing:
        try:
            from slab_utils import trimesh_to_freecad
            wing_fc_mesh = trimesh_to_freecad(wing_mesh)
            if wing_fc_mesh:
                from viz_utils import show_mesh
                show_mesh(wing_fc_mesh, doc, "WingMesh",
                          color=(0.5, 0.5, 0.5), transparency=50)
                doc.recompute()
                print("Wing mesh visualised.")
        except Exception as e:
            print(f"Wing mesh visualisation error: {e}")

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
        obj_name='WingR1_msv_orient001_solid',
        rib_spacing=20.0,
        rib_angle=30.0,
        grid_orientation=0.0,
        primary_dir=FreeCAD.Vector(0, 0, 1),
        bridge_height=0.4,
        hole_margin=0.5,
        thickness=0.4,
        input_step_path="",
        vis_rib_centre_surfaces=False,
        vis_rib_centre_surfaces_clip=False,
        vis_rib_segments=True,
        vis_centre_lines=False,
        vis_bridge=True,
        vis_hole=True,
        vis_wing=True,
        z_step=0.2
    )

    params.plane_normal = pdef['normal']
    params.plane_axis_u = pdef['axis_u']
    params.plane_axis_v = pdef['axis_v']

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