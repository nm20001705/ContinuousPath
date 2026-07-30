# bridge_utils.py

import FreeCAD
import Part
import Mesh
import MeshPart
import numpy as np
import trimesh
from shapely.geometry import Polygon as ShapelyPolygon, LineString
import trimesh.path.polygons
from shapely.geometry import Polygon as ShapelyPolygon, LineString, Point
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad
from shapely.ops import linemerge, polygonize
from scipy.spatial import ConvexHull

def create_bridges_analytical(wing_mesh, rib_segments, primary_dir_np,
                              z_step=0.2, bridge_height=0.5,
                              z_vals=None, slices=None,
                              doc=None, vis=False):
    prim = primary_dir_np / np.linalg.norm(primary_dir_np)

    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)

    all_vertices = []
    all_faces = []
    vert_offset = 0

    # ----- Debug counters -----
    total_slices_visited = 0
    skipped_chord_empty = 0
    skipped_degenerate = 0
    skipped_chord_len = 0
    skipped_invalid_poly = 0
    skipped_midpoint = 0
    skipped_seg_width = 0
    kept = 0

    for seg in rib_segments:
        p0 = seg['p0']
        p1 = seg['p1']
        dir_rib = seg['dir']
        slab_normal = seg['slab_normal']

        d0 = np.dot(p0, prim)
        d1 = np.dot(p1, prim)
        d_min = min(d0, d1)
        d_max = max(d0, d1)

        line_dir_3d = np.cross(slab_normal, prim)
        norm_l = np.linalg.norm(line_dir_3d)
        if norm_l < 1e-8:
            continue
        line_dir_3d /= norm_l
        bridge_dir_2d = np.array([np.dot(line_dir_3d, u_ax), np.dot(line_dir_3d, v_ax)])
        bridge_dir_2d /= np.linalg.norm(bridge_dir_2d)

        ribbon_pts = []

        start_idx = np.searchsorted(z_vals, d_min)
        end_idx = np.searchsorted(z_vals, d_max, side='right')

        for idx in range(start_idx, end_idx):
            d = z_vals[idx]
            wing_poly, to_3d_mat = slices[idx]
            total_slices_visited += 1

            denom = np.dot(prim, dir_rib)
            if abs(denom) < 1e-8:
                continue
            t_plane = np.dot(prim, d * prim - p0) / denom
            P_3d = p0 + t_plane * dir_rib

            P_uv = np.array([np.dot(P_3d - d * prim, u_ax),
                             np.dot(P_3d - d * prim, v_ax)])

            wing_uv = ShapelyPolygon(np.array(wing_poly.exterior.coords))

            extent = 1e5
            line = LineString([P_uv - extent * bridge_dir_2d, P_uv + extent * bridge_dir_2d])
            try:
                chord = wing_uv.intersection(line)
            except Exception:
                continue
            if chord.is_empty:
                skipped_chord_empty += 1
                continue

            if chord.geom_type == 'MultiLineString':
                best_seg = None
                min_dist = np.inf
                for sub in chord.geoms:
                    dist = sub.distance(Point(P_uv))
                    if dist < min_dist:
                        min_dist = dist
                        best_seg = sub
                chord = best_seg
            if chord is None or chord.is_empty:
                skipped_chord_empty += 1
                continue

            coords = np.array(chord.coords)
            if len(coords) < 2:
                skipped_degenerate += 1
                continue

            chord_len = np.linalg.norm(coords[-1] - coords[0])
            if chord_len < 0.1:
                skipped_chord_len += 1
                continue

            if not wing_uv.is_valid:
                skipped_invalid_poly += 1
                continue

            midpoint_uv = (coords[0] + coords[-1]) / 2.0
            if not wing_uv.buffer(1e-6).contains(Point(midpoint_uv)):
                continue

            t_vals = np.dot(coords, bridge_dir_2d)
            t_min, t_max = t_vals.min(), t_vals.max()
            t_center = np.dot(midpoint_uv, bridge_dir_2d)

            half_w = bridge_height
            t_start = max(t_center - half_w, t_min)
            t_end   = min(t_center + half_w, t_max)

            if t_end - t_start < 1e-6:
                skipped_seg_width += 1
                continue

            seg_start_uv = midpoint_uv + (t_start - t_center) * bridge_dir_2d
            seg_end_uv   = midpoint_uv + (t_end   - t_center) * bridge_dir_2d

            def uv_to_3d(uv):
                pt = np.array([uv[0], uv[1], 0.0, 1.0])
                pt3d = to_3d_mat @ pt
                return pt3d[:3]

            pt_left  = uv_to_3d(seg_start_uv)
            pt_right = uv_to_3d(seg_end_uv)
            ribbon_pts.append((d, pt_left, pt_right))
            kept += 1

        if len(ribbon_pts) >= 2:
            ribbon_pts.sort(key=lambda x: x[0])
            for k in range(len(ribbon_pts) - 1):
                _, L1, R1 = ribbon_pts[k]
                _, L2, R2 = ribbon_pts[k+1]
                v0 = vert_offset
                v1 = vert_offset + 1
                v2 = vert_offset + 2
                v3 = vert_offset + 3
                all_vertices.extend([L1, R1, R2, L2])
                all_faces.append([v0, v1, v2])
                all_faces.append([v0, v2, v3])
                vert_offset += 4

    print(f"DEBUG: slices visited: {total_slices_visited}")
    print(f"  chord empty: {skipped_chord_empty}")
    print(f"  degenerate (len<2): {skipped_degenerate}")
    print(f"  chord_len < 0.1: {skipped_chord_len}")
    print(f"  invalid polygon: {skipped_invalid_poly}")
    print(f"  midpoint outside: {skipped_midpoint}")
    print(f"  segment width zero: {skipped_seg_width}")
    print(f"  kept: {kept}")

    if not all_vertices:
        print("No bridge geometry generated.")
        return None

    # ... rest unchanged ...
