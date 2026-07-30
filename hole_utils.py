# hole_utils.py

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

def create_holes(wing_mesh, rib_segments, primary_dir_np,
                 z_step=0.2, point_condition=None, hole_margin=0.0,
                 doc=None, vis=False):
    """
    Build hole ribbons inside each rib segment, with a safety margin from
    the wing surface and from the top/bottom of the segment.

    Parameters
    ----------
    hole_margin : float
        Distance (in model units) to keep away from the wing surface (on both
        sides of the chord) and from the top/bottom of the segment's Z extent.
    """
    if point_condition is None:
        point_condition = lambda x: 0.0

    prim = primary_dir_np / np.linalg.norm(primary_dir_np)

    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)

    def wing_poly_at_plane(plane_origin):
        # ... identical to previous version ...
        try:
            result = trimesh.intersections.mesh_plane(
                wing_mesh, plane_normal=prim, plane_origin=plane_origin
            )
        except Exception:
            return None
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
        to_3d = np.eye(4)
        to_3d[:3, 0] = u_ax
        to_3d[:3, 1] = v_ax
        to_3d[:3, 3] = origin
        return poly, to_3d

    all_vertices = []
    all_faces = []
    vert_offset = 0
    num_holes = 0
    total_segments = len(rib_segments)

    for seg in rib_segments:
        p0 = seg['p0']
        p1 = seg['p1']
        dir_rib = seg['dir']
        slab_normal = seg['slab_normal']

        d0 = np.dot(p0, prim)
        d1 = np.dot(p1, prim)
        d_min = min(d0, d1)
        d_max = max(d0, d1)

        # Apply Z margin
        d_start = d_min + hole_margin
        d_end   = d_max - hole_margin
        if d_start >= d_end:
            continue
        segment_z_span = d_end - d_start   # used for normalised position

        line_dir_3d = np.cross(slab_normal, prim)
        norm_l = np.linalg.norm(line_dir_3d)
        if norm_l < 1e-8:
            continue
        line_dir_3d /= norm_l
        bridge_dir_2d = np.array([np.dot(line_dir_3d, u_ax), np.dot(line_dir_3d, v_ax)])
        bridge_dir_2d /= np.linalg.norm(bridge_dir_2d)

        ribbon_pts = []
        d = d_start
        while d <= d_end + 1e-9:
            plane_origin = d * prim
            res = wing_poly_at_plane(plane_origin)
            if res is None:
                d += z_step
                continue
            wing_poly, to_3d_mat = res

            denom = np.dot(prim, dir_rib)
            if abs(denom) < 1e-8:
                d += z_step
                continue
            t_plane = np.dot(prim, plane_origin - p0) / denom
            P_3d = p0 + t_plane * dir_rib

            P_uv = np.array([np.dot(P_3d - plane_origin, u_ax),
                             np.dot(P_3d - plane_origin, v_ax)])

            wing_uv = ShapelyPolygon(np.array(wing_poly.exterior.coords))
            if not wing_uv.contains(Point(P_uv)):
                d += z_step
                continue

            extent = 1e5
            line = LineString([P_uv - extent * bridge_dir_2d, P_uv + extent * bridge_dir_2d])
            try:
                chord = wing_uv.intersection(line)
            except Exception:
                d += z_step
                continue
            if chord.is_empty:
                d += z_step
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
                d += z_step
                continue

            coords = np.array(chord.coords)
            t_vals = np.dot(coords, bridge_dir_2d)
            t_min, t_max = t_vals.min(), t_vals.max()

            # Shrink chord by hole_margin on each side
            t_min_safe = t_min + hole_margin
            t_max_safe = t_max - hole_margin
            if t_min_safe >= t_max_safe:
                d += z_step
                continue

            t_center = np.dot(P_uv, bridge_dir_2d)
            chord_half_effective = (t_max_safe - t_min_safe) / 2.0

            if segment_z_span > 1e-8:
                x = (d - d_start) / segment_z_span   # 0..1
            else:
                x = 0.0
            factor = point_condition(x)
            half_w = factor * chord_half_effective

            t_start = max(t_center - half_w, t_min_safe)
            t_end   = min(t_center + half_w, t_max_safe)

            seg_start_uv = P_uv + (t_start - t_center) * bridge_dir_2d
            seg_end_uv   = P_uv + (t_end   - t_center) * bridge_dir_2d

            def uv_to_3d(uv):
                pt = np.array([uv[0], uv[1], 0.0, 1.0])
                pt3d = to_3d_mat @ pt
                return pt3d[:3]

            pt_left  = uv_to_3d(seg_start_uv)
            pt_right = uv_to_3d(seg_end_uv)
            ribbon_pts.append((d, pt_left, pt_right))
            d += z_step

        if len(ribbon_pts) >= 2:
            num_holes += 1
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

    if not all_vertices:
        print("No hole geometry generated.")
        return None

    verts_arr = np.array(all_vertices)
    faces_arr = np.array(all_faces)
    hole_mesh = trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False)
    hole_mesh.merge_vertices()
    hole_mesh.fix_normals()

    print(f"Hole mesh: {len(hole_mesh.vertices)} verts, {len(hole_mesh.faces)} faces")
    print(f"Holes created: {num_holes} out of {total_segments} rib segments")

    if vis and doc:
        try:
            fc_mesh = trimesh_to_freecad(hole_mesh)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "Holes",
                          color=(0.1, 0.1, 0.1), transparency=80)
                doc.recompute()
                print("Holes visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return hole_mesh
