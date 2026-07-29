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
from shapely.geometry import Polygon as ShapelyPolygon, LineString
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad

def create_bridges_analytical(wing_mesh, all_lines_np, plane_normal_np,
                              z_step=0.2, bridge_height=0.5,
                              doc=None, vis=False):
    """
    Create bridge surfaces aligned with the rib direction inside each rib segment.

    For each cell (between two rib crossings), a narrow rectangular strip is
    placed centred on the cell's midline, oriented along the rib direction,
    and intersected with the wing cross‑section polygon.
    The result is an elongated ribbon that follows the rib.

    Parameters
    ----------
    wing_mesh : trimesh.Trimesh
    all_lines_np : list of (start, end) numpy arrays
    plane_normal_np : (3,) np.array
    z_step : float          (unused in this version, kept for compatibility)
    bridge_height : float   half‑width of the bridge (perpendicular to rib direction)
    doc, vis : as before.
    """
    # ---- Precompute line data ----
    lines_data = []
    for (start, end) in all_lines_np:
        d = end - start
        norm = np.linalg.norm(d)
        if norm < 1e-8:
            lines_data.append(None)
            continue
        d /= norm
        lines_data.append({'point': start, 'dir': d})

    def intersect_lines(p1, d1, p2, d2):
        A = np.column_stack((d1, -d2))
        b = p2 - p1
        try:
            st, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)
        except np.linalg.LinAlgError:
            return None
        if rank < 2:
            return None
        s, t = st[0], st[1]
        p = p1 + s * d1
        if np.linalg.norm(p2 + t * d2 - p) > 1e-6:
            return None
        return p

    all_bridge_meshes = []

    for i, (line_i, ld_i) in enumerate(zip(all_lines_np, lines_data)):
        if ld_i is None:
            continue

        p_i = ld_i['point']
        d_i = ld_i['dir']          # unit rib direction (3D)

        # ---- Rib plane normal ----
        N_i = np.cross(d_i, plane_normal_np)
        if np.linalg.norm(N_i) < 1e-8:
            continue
        N_i /= np.linalg.norm(N_i)

        # ---- 1. Wing cross‑section polygon + local 2D transform ----
        try:
            section = wing_mesh.section(plane_origin=p_i, plane_normal=N_i)
        except Exception:
            continue
        if section is None:
            continue
        try:
            planar, to_3d = section.to_planar()
        except Exception:
            continue
        polygons = planar.polygons_full
        if not polygons:
            continue
        wing_poly = polygons[0]          # outer boundary
        if wing_poly.is_empty or wing_poly.area < 1e-8:
            continue

        # Local 2D basis from the transform
        to_3d_mat = np.array(to_3d)      # (4,4)
        u_ax = to_3d_mat[:3, 0]          # local X axis (3D)
        v_ax = to_3d_mat[:3, 1]          # local Y axis (3D)
        origin_3d = to_3d_mat[:3, 3]

        # ---- 2. Project wing polygon and rib direction into local 2D ----
        xy = np.array(wing_poly.exterior.coords)        # (N,2)
        ones = np.ones((len(xy), 1))
        xy_h = np.hstack([xy, np.zeros((len(xy), 1)), ones])
        pts_3d = (to_3d_mat @ xy_h.T).T[:, :3]

        # Map polygon to 2D
        u_wing = np.dot(pts_3d - origin_3d, u_ax)
        v_wing = np.dot(pts_3d - origin_3d, v_ax)
        wing_uv = ShapelyPolygon(np.column_stack([u_wing, v_wing]))
        if wing_uv.is_empty or wing_uv.area < 1e-8:
            continue

        # Project the rib direction d_i into 2D
        rib_dir_2d = np.array([np.dot(d_i, u_ax), np.dot(d_i, v_ax)])
        norm_2d = np.linalg.norm(rib_dir_2d)
        if norm_2d < 1e-8:
            continue
        rib_dir_2d /= norm_2d
        # Perpendicular direction (for bridge half‑width)
        perp_2d = np.array([-rib_dir_2d[1], rib_dir_2d[0]])

        # ---- 3. Find crossing u‑coordinates (in local 2D, using rib_dir_2d) ----
        # We need the position along rib_dir_2d for crossing lines.
        crossing_t = []  # parameter along rib_dir_2d
        for j, ld_j in enumerate(lines_data):
            if j == i or ld_j is None:
                continue
            if abs(np.dot(d_i, ld_j['dir'])) > 0.9999:
                continue
            pt = intersect_lines(p_i, d_i, ld_j['point'], ld_j['dir'])
            if pt is None:
                continue
            # 2D coordinates of the intersection point
            uv_pt = np.array([np.dot(pt - origin_3d, u_ax),
                              np.dot(pt - origin_3d, v_ax)])
            # Project onto rib_dir_2d to get the t parameter
            t = np.dot(uv_pt, rib_dir_2d)
            crossing_t.append(t)

        # Get the t-range of the wing polygon
        t_wing = np.dot(np.column_stack([u_wing, v_wing]), rib_dir_2d)
        t_min = t_wing.min() - 10.0
        t_max = t_wing.max() + 10.0
        crossing_t.extend([t_min, t_max])
        crossing_t = sorted(set(crossing_t))

        # ---- 4. For each interval along rib_dir, build a bridge strip ----
        for k in range(len(crossing_t)-1):
            t0 = crossing_t[k]
            t1 = crossing_t[k+1]
            if t1 - t0 < 1e-8:
                continue

            # Cell polygon = wing_uv ∩ strip bounded by lines t = t0 and t = t1
            # Define the two bounding lines
            # Line perpendicular to rib_dir_2d, passing through points with t=t0,t1
            # point on line: P = t0 * rib_dir_2d + s * perp_2d
            # We can create a huge rectangle that covers the cell extent in perp direction.
            # Get perpendicular extent of wing_uv
            perp_wing = np.dot(np.column_stack([u_wing, v_wing]), perp_2d)
            p_min = perp_wing.min() - 10.0
            p_max = perp_wing.max() + 10.0

            # Rectangle aligned with rib_dir_2d, spanning [t0,t1] along rib and [p_min,p_max] along perp
            # Four corners in order: (t0,p_min), (t1,p_min), (t1,p_max), (t0,p_max)
            corners = [
                (t0, p_min), (t1, p_min), (t1, p_max), (t0, p_max)
            ]
            # Transform to (u,v) coordinates: (u,v) = t*rib_dir_2d + p*perp_2d
            rect_pts = []
            for (t_val, p_val) in corners:
                pt = t_val * rib_dir_2d + p_val * perp_2d
                rect_pts.append(pt)
            rect = ShapelyPolygon(rect_pts)

            try:
                cell = wing_uv.intersection(rect)
            except Exception:
                continue
            if cell.is_empty:
                continue
            if cell.geom_type == 'Polygon':
                cells = [cell]
            elif cell.geom_type == 'MultiPolygon':
                cells = list(cell.geoms)
            else:
                continue

            for c in cells:
                if c.is_empty or c.area < 1e-8:
                    continue
                # Build the bridge strip: narrow rectangle centred on the midline t_mid = (t0+t1)/2
                t_mid = (t0 + t1) / 2.0
                # Half‑width along rib_dir_2d is bridge_height (since bridge_height is half-length along the rib direction)
                # But the cell width is t1 - t0; we want a strip narrower than the cell.
                # Use bridge_height as the half‑width along the rib direction.
                half_w = min(bridge_height, (t1 - t0) * 0.4)  # ensure it doesn't exceed cell
                t_left = t_mid - half_w
                t_right = t_mid + half_w

                # Rectangle for the bridge strip: aligned with rib_dir_2d, from t_left to t_right
                strip_corners = [
                    (t_left, p_min), (t_right, p_min), (t_right, p_max), (t_left, p_max)
                ]
                strip_pts = []
                for (t_val, p_val) in strip_corners:
                    pt = t_val * rib_dir_2d + p_val * perp_2d
                    strip_pts.append(pt)
                strip_rect = ShapelyPolygon(strip_pts)

                # Intersect the strip with the cell
                try:
                    bridge_poly = c.intersection(strip_rect)
                except Exception:
                    continue
                if bridge_poly.is_empty:
                    continue

                if bridge_poly.geom_type == 'Polygon':
                    bridges = [bridge_poly]
                elif bridge_poly.geom_type == 'MultiPolygon':
                    bridges = list(bridge_poly.geoms)
                else:
                    continue

                for bp in bridges:
                    if bp.is_empty or bp.area < 1e-8:
                        continue
                    try:
                        verts2d, faces = trimesh.creation.triangulate_polygon(bp)
                    except Exception:
                        continue
                    if len(verts2d) == 0 or len(faces) == 0:
                        continue

                    verts_h = np.column_stack([verts2d, np.zeros(len(verts2d)), np.ones(len(verts2d))])
                    verts3d = (to_3d_mat @ verts_h.T).T[:, :3]
                    mesh = trimesh.Trimesh(vertices=verts3d, faces=faces, process=False)
                    mesh.fix_normals()
                    all_bridge_meshes.append(mesh)

    print(f"Created {len(all_bridge_meshes)} bridge triangles.")

    if not all_bridge_meshes:
        return None

    bridges_mesh = trimesh.util.concatenate(all_bridge_meshes)

    if vis and doc and bridges_mesh:
        try:
            fc_mesh = trimesh_to_freecad(bridges_mesh)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "Bridges",
                          color=(0.9, 0.7, 0.1), transparency=20)
                doc.recompute()
                print("Bridges visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return bridges_mesh
