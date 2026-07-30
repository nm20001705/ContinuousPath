# bridge_utils.py

import FreeCAD
import Part
import Mesh
import MeshPart
import math
import numpy as np
import trimesh
import shapely.geometry as sg
from shapely.geometry import Polygon as ShapelyPolygon, LineString
import trimesh.path.polygons
from shapely.geometry import Polygon as ShapelyPolygon, LineString, Point
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad
from shapely.ops import linemerge, polygonize
from scipy.spatial import ConvexHull

def create_bridges_analytical(wing_mesh, all_lines_np, primary_dir_np,
                              z_step=0.2, bridge_height=0.5,
                              doc=None, vis=False):
    """
    Build bridge ribbons inside each slab by slicing the wing perpendicular to
    primary_dir_np, finding the rib chord at every slice, and connecting narrow
    central segments into a mesh.
    """
    prim = primary_dir_np / np.linalg.norm(primary_dir_np)

    # ---- Precompute rib line data ----
    lines_data = []
    for (start, end) in all_lines_np:
        d = end - start
        norm = np.linalg.norm(d)
        if norm < 1e-8:
            lines_data.append(None)
            continue
        d /= norm
        lines_data.append({'point': start, 'dir': d})

    print(f"Total rib lines: {len(all_lines_np)}, valid: {sum(ld is not None for ld in lines_data)}")

    # ---- Slice range from wing bounding box ----
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
    proj = np.dot(corners, prim)
    d_min = proj.min() - 1e-3
    d_max = proj.max() + 1e-3
    print(f"Slice range: d_min={d_min:.3f}, d_max={d_max:.3f}, z_step={z_step}")

    # ---- Fixed basis vectors for all slice planes ----
    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)

    def wing_poly_at_plane(plane_origin):
        """
        Returns (2D polygon, 4x4 transform) of the wing cross‑section at the given plane.
        Uses trimesh.intersections.mesh_plane + Shapely polygonization.
        """
        # ---------- FIX: safe unpacking for all trimesh versions ----------
        try:
            result = trimesh.intersections.mesh_plane(
                wing_mesh, plane_normal=prim, plane_origin=plane_origin
            )
        except Exception as e:
            print(f"  [mesh_plane exception] origin={plane_origin} error={e}")
            return None

        # result can be (lines,) or (lines, face_index) or (lines, face_index, valid)
        if isinstance(result, tuple):
            lines = result[0]
        else:
            lines = result

        if lines is None or len(lines) == 0:
            return None

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
            return None

        merged = linemerge(segments_2d)
        if merged.is_empty:
            return None

        polys = list(polygonize(merged))
        if not polys:
            # Fallback to convex hull of all endpoints
            all_pts = []
            for line in segments_2d:
                all_pts.extend(line.coords)
            pts = np.array(all_pts)
            if len(pts) < 3:
                return None
            try:
                hull = ConvexHull(pts)
                poly = ShapelyPolygon(pts[hull.vertices])
            except Exception:
                return None
        else:
            poly = max(polys, key=lambda p: p.area)

        if poly.is_empty or poly.area < 1e-8:
            return None

        # Build 4x4 transform (u_ax, v_ax, origin)
        to_3d = np.eye(4)
        to_3d[:3, 0] = u_ax
        to_3d[:3, 1] = v_ax
        to_3d[:3, 3] = origin

        return poly, to_3d

    # ---- Main loop over ribs ----
    all_vertices = []
    all_faces = []
    vert_offset = 0
    total_slices, valid_slices = 0, 0

    for i, ld_i in enumerate(lines_data):
        if ld_i is None:
            continue

        p0 = ld_i['point']
        dir_rib = ld_i['dir']
        print(f"\nRib {i}: p0={p0}, dir={dir_rib}")

        ribbon_pts = []
        d = d_min
        while d <= d_max + 1e-9:
            # plane position is at coordinate d along the primary direction
            plane_origin = d * prim
            total_slices += 1

            res = wing_poly_at_plane(plane_origin)
            if res is None:
                d += z_step
                continue
            wing_poly, to_3d_mat = res

            # Convert polygon to 2D local coordinates (u, v)
            xy = np.array(wing_poly.exterior.coords)
            wing_uv = ShapelyPolygon(xy)

            # Rib direction projected to 2D
            rib_dir_2d = np.array([np.dot(dir_rib, u_ax), np.dot(dir_rib, v_ax)])
            norm_2d = np.linalg.norm(rib_dir_2d)
            if norm_2d < 1e-8:
                d += z_step
                continue
            rib_dir_2d /= norm_2d

            # Intersection of the 3D rib line with the slice plane
            denom = np.dot(prim, dir_rib)
            if abs(denom) < 1e-8:
                d += z_step
                continue
            t_plane = np.dot(prim, plane_origin - p0) / denom
            P_3d = p0 + t_plane * dir_rib

            # 2D coordinates of the rib intersection point
            P_uv = np.array([np.dot(P_3d - plane_origin, u_ax),
                             np.dot(P_3d - plane_origin, v_ax)])

            # Find chord through P_uv along rib_dir_2d
            extent = 1e5
            line = LineString([P_uv - extent * rib_dir_2d, P_uv + extent * rib_dir_2d])
            try:
                chord = wing_uv.intersection(line)
            except Exception:
                d += z_step
                continue
            if chord.is_empty:
                d += z_step
                continue

            # If multiple segments, pick the one containing P
            if chord.geom_type == 'MultiLineString':
                best_seg = None
                min_dist = np.inf
                for seg in chord.geoms:
                    dist = seg.distance(Point(P_uv))
                    if dist < min_dist:
                        min_dist = dist
                        best_seg = seg
                chord = best_seg
            if chord is None or chord.is_empty:
                d += z_step
                continue

            # Chord endpoints and the position of P along the direction
            coords = np.array(chord.coords)
            t_vals = np.dot(coords, rib_dir_2d)
            t_min, t_max = t_vals.min(), t_vals.max()
            t_center = np.dot(P_uv, rib_dir_2d)

            # Build bridge segment centered on the rib intersection, clipped to chord
            half_w = bridge_height
            t_left  = max(t_center - half_w, t_min)
            t_right = min(t_center + half_w, t_max)

            seg_start_uv = t_left  * rib_dir_2d
            seg_end_uv   = t_right * rib_dir_2d

            # 2D → 3D
            def uv_to_3d(uv):
                pt = np.array([uv[0], uv[1], 0.0, 1.0])
                pt3d = to_3d_mat @ pt
                return pt3d[:3]

            pt_left  = uv_to_3d(seg_start_uv)
            pt_right = uv_to_3d(seg_end_uv)

            ribbon_pts.append((d, pt_left, pt_right))
            valid_slices += 1
            d += z_step

        print(f"  Rib {i}: collected {len(ribbon_pts)} ribbon points")

        # Build quad mesh for this rib
        if len(ribbon_pts) >= 2:
            ribbon_pts.sort(key=lambda x: x[0])  # sort by d
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

    print(f"\nTotal slices: {total_slices}, valid: {valid_slices}")
    print(f"Total vertices: {len(all_vertices)}")

    if not all_vertices:
        print("No bridge geometry generated.")
        return None

    verts_arr = np.array(all_vertices)
    faces_arr = np.array(all_faces)
    bridge_mesh = trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False)
    bridge_mesh.merge_vertices()
    bridge_mesh.fix_normals()

    print(f"Bridge mesh: {len(bridge_mesh.vertices)} verts, {len(bridge_mesh.faces)} faces")

    if vis and doc:
        try:
            fc_mesh = trimesh_to_freecad(bridge_mesh)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "Bridges",
                          color=(0.9, 0.7, 0.1), transparency=20)
                doc.recompute()
                print("Bridges visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return bridge_mesh
