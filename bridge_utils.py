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

# bridge_utils.py – corrected bridge direction and centering

import trimesh
import numpy as np
from shapely.geometry import LineString, Polygon as ShapelyPolygon, Point
from shapely.ops import linemerge, polygonize
from scipy.spatial import ConvexHull
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad

def create_bridges_analytical(wing_mesh, all_lines_np, primary_dir_np,
                              z_step=0.2, bridge_height=0.5,
                              construction_plane_normal=None,
                              doc=None, vis=False):
    """
    Build bridge ribbons inside each slab by slicing the wing perpendicular to
    primary_dir_np, using the true slab‑slice intersection line, and centering
    narrow segments on the rib intersection point.
    """
    prim = primary_dir_np / np.linalg.norm(primary_dir_np)

    if construction_plane_normal is None:
        raise ValueError("construction_plane_normal must be provided")
    con_normal = construction_plane_normal / np.linalg.norm(construction_plane_normal)

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
    d_min = proj.min() - z_step * 0.5
    d_max = proj.max() + z_step * 0.5

    # ---- Fixed basis for slice plane (perpendicular to prim) ----
    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)

    def wing_poly_at_plane(plane_origin):
        """Return (2D polygon, 4x4 transform) or None."""
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

    # ---- Main loop over ribs ----
    all_vertices = []
    all_faces = []
    vert_offset = 0

    for i, ld_i in enumerate(lines_data):
        if ld_i is None:
            continue

        p0 = ld_i['point']
        dir_rib = ld_i['dir']

        # ---- Slab plane normal: perpendicular to rib direction and construction normal ----
        slab_normal = np.cross(dir_rib, con_normal)
        norm_s = np.linalg.norm(slab_normal)
        if norm_s < 1e-8:
            continue
        slab_normal /= norm_s

        # ---- Direction of slab‑slice intersection line (in 3D) ----
        line_dir_3d = np.cross(slab_normal, prim)
        norm_l = np.linalg.norm(line_dir_3d)
        if norm_l < 1e-8:
            continue
        line_dir_3d /= norm_l

        # Project into the slice plane’s 2D coordinates
        bridge_dir_2d = np.array([np.dot(line_dir_3d, u_ax), np.dot(line_dir_3d, v_ax)])
        bridge_dir_2d /= np.linalg.norm(bridge_dir_2d)

        ribbon_pts = []
        d = d_min
        while d <= d_max + 1e-9:
            plane_origin = d * prim
            res = wing_poly_at_plane(plane_origin)
            if res is None:
                d += z_step
                continue
            wing_poly, to_3d_mat = res

            # ---- Intersection of rib line with this slice plane ----
            denom = np.dot(prim, dir_rib)
            if abs(denom) < 1e-8:
                d += z_step
                continue
            t_plane = np.dot(prim, plane_origin - p0) / denom
            P_3d = p0 + t_plane * dir_rib

            # 2D coordinates of the rib intersection point
            P_uv = np.array([np.dot(P_3d - plane_origin, u_ax),
                             np.dot(P_3d - plane_origin, v_ax)])

            # Only keep this slice if the rib point lies inside the wing cross‑section
            wing_uv = ShapelyPolygon(np.array(wing_poly.exterior.coords))
            if not wing_uv.contains(Point(P_uv)):
                d += z_step
                continue

            # Find the chord along bridge_dir_2d through P_uv
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
                for seg in chord.geoms:
                    dist = seg.distance(Point(P_uv))
                    if dist < min_dist:
                        min_dist = dist
                        best_seg = seg
                chord = best_seg
            if chord is None or chord.is_empty:
                d += z_step
                continue

            coords = np.array(chord.coords)
            # Parameter range along bridge_dir_2d
            t_vals = np.dot(coords, bridge_dir_2d)
            t_min, t_max = t_vals.min(), t_vals.max()
            t_center = np.dot(P_uv, bridge_dir_2d)

            # Bridge segment centered on P_uv, clipped to chord
            half_w = bridge_height
            t_start = max(t_center - half_w, t_min)
            t_end   = min(t_center + half_w, t_max)

            seg_start_uv = P_uv + (t_start - t_center) * bridge_dir_2d
            seg_end_uv   = P_uv + (t_end   - t_center) * bridge_dir_2d

            # 2D → 3D
            def uv_to_3d(uv):
                pt = np.array([uv[0], uv[1], 0.0, 1.0])
                pt3d = to_3d_mat @ pt
                return pt3d[:3]

            pt_left  = uv_to_3d(seg_start_uv)
            pt_right = uv_to_3d(seg_end_uv)

            ribbon_pts.append((d, pt_left, pt_right))
            d += z_step

        # Build quad mesh for this rib
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
