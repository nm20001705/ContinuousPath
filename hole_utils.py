# hole_utils.py – analytical version using rib_segment meshes directly
#
# Instead of intersecting the FULL wing slice polygon (which can contain
# internal holes -> polygon-with-holes headaches), we slice each rib
# segment mesh directly. A rib segment is already a solid chunk of the
# wing (the result of a boolean intersection), so any horizontal slice of
# it is guaranteed to be material only.
#
# VOID HANDLING:
# If a rib segment slice produces multiple disjoint line pieces (e.g. the
# cell straddles a servo-hole void), we cluster the raw mesh_plane edges
# into connected pieces (by shared endpoints) and pick whichever piece is
# CLOSEST to the rib centerline's projection (t_rib) -- distance 0 if
# t_rib genuinely falls inside it. hole_margin is then applied to that
# piece's own boundaries -- which may be the outer wing skin on one side
# and the void wall on the other. This makes the hole shrink away from
# the void by hole_margin automatically, with no special-casing needed:
# a void boundary is just another edge.
#
# Nearest-piece (rather than strict containment) matters most near thin
# sections (e.g. the trailing edge): t_rib comes from an idealized
# infinite centerline, while piece bounds come from the triangulated
# mesh, and once the solid strip gets thinner than that numerical noise,
# strict containment spuriously rejects perfectly good slices. A slice is
# only skipped now if mesh_plane returns literally nothing to cluster.

import numpy as np
import trimesh
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad


def _cluster_segments(pts_a, pts_b, tol=5):
    """
    Group segment indices into connected components based on shared
    endpoints. pts_a/pts_b are (N,2) arrays of uv endpoints for N segments.
    Returns a list of index-arrays, one per connected piece.
    """
    n = len(pts_a)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    key_map = {}
    for idx in range(n):
        for endpoint in (pts_a[idx], pts_b[idx]):
            key = (round(float(endpoint[0]), tol), round(float(endpoint[1]), tol))
            if key in key_map:
                union(idx, key_map[key])
            else:
                key_map[key] = idx

    groups = {}
    for idx in range(n):
        r = find(idx)
        groups.setdefault(r, []).append(idx)

    return list(groups.values())


def create_holes_analytical(rib_segment_meshes, bridge_segments, primary_dir_np,
                             z_vals, point_condition=None, hole_margin=0.0,
                             doc=None, vis=False):
    """
    rib_segment_meshes : list[trimesh.Trimesh]
        The solid rib cell pieces from build_rib_segments_analytical().
    bridge_segments : list[dict]
        Same length / same order as rib_segment_meshes. Each dict has
        'p0', 'p1', 'dir', 'slab_normal' (as built in main.py from
        segment_bounds). The i-th bridge_segment corresponds to the i-th
        rib_segment_meshes entry because both come from the same
        segment_bounds list, in the same order.
    z_vals : np.ndarray
        Global precomputed Z-slice positions (used only to keep slice
        spacing consistent with the rest of the pipeline).
    """
    if point_condition is None:
        point_condition = lambda x: 1.0

    if len(rib_segment_meshes) != len(bridge_segments):
        raise ValueError(
            f"rib_segment_meshes ({len(rib_segment_meshes)}) and "
            f"bridge_segments ({len(bridge_segments)}) must be the same "
            f"length and correspond index-for-index."
        )

    prim = primary_dir_np / np.linalg.norm(primary_dir_np)

    # Same global basis used everywhere else in the pipeline
    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)

    all_vertices = []
    all_faces = []
    vert_offset = 0

    for seg_mesh, seg in zip(rib_segment_meshes, bridge_segments):
        if seg_mesh is None or len(seg_mesh.vertices) == 0:
            continue

        p0 = seg['p0']
        dir_rib = seg['dir']
        slab_normal = seg['slab_normal']

        # Bridge direction (perpendicular to rib centerline, in-slice)
        line_dir_3d = np.cross(slab_normal, prim)
        norm_l = np.linalg.norm(line_dir_3d)
        if norm_l < 1e-8:
            continue
        line_dir_3d /= norm_l
        bridge_dir_2d = np.array([np.dot(line_dir_3d, u_ax),
                                   np.dot(line_dir_3d, v_ax)])
        bd_norm = np.linalg.norm(bridge_dir_2d)
        if bd_norm < 1e-8:
            continue
        bridge_dir_2d /= bd_norm

        def uv_of(pt3d, d):
            delta = pt3d - d * prim
            return np.array([np.dot(delta, u_ax), np.dot(delta, v_ax)])

        def t_of(uv):
            return np.dot(uv, bridge_dir_2d)

        # Z-range of this specific rib segment mesh (projected onto prim)
        proj = seg_mesh.vertices @ prim
        d_min = proj.min()
        d_max = proj.max()
        if d_max - d_min < 1e-8:
            continue

        start_idx = np.searchsorted(z_vals, d_min)
        end_idx = np.searchsorted(z_vals, d_max, side='right')

        ribbons = []  # (d, left_pt3d, right_pt3d)

        for idx in range(start_idx, end_idx):
            d = z_vals[idx]
            plane_origin = d * prim

            try:
                result = trimesh.intersections.mesh_plane(
                    seg_mesh, plane_normal=prim, plane_origin=plane_origin
                )
            except Exception:
                continue

            lines = result[0] if isinstance(result, tuple) else result
            if lines is None or len(lines) == 0:
                continue

            lines_arr = np.asarray(lines)
            if len(lines_arr) < 1:
                continue

            # uv endpoints of every raw segment returned by mesh_plane
            pts_a_uv = np.array([uv_of(p, d) for p in lines_arr[:, 0, :]])
            pts_b_uv = np.array([uv_of(p, d) for p in lines_arr[:, 1, :]])

            # Rib centerline position on this plane, used both to pick
            # the correct piece and as the reference origin for
            # reconstructing points on the straight bridge line.
            denom = np.dot(prim, dir_rib)
            if abs(denom) < 1e-8:
                continue
            t_plane = np.dot(prim, plane_origin - p0) / denom
            P_3d = p0 + t_plane * dir_rib
            P_uv = uv_of(P_3d, d)
            t_rib = t_of(P_uv)

            # Cluster raw segments into connected pieces (handles voids:
            # a piece boundary can be the outer skin OR a void wall).
            groups = _cluster_segments(pts_a_uv, pts_b_uv)

            # Pick the piece CLOSEST to t_rib rather than requiring strict
            # containment. Near the trailing edge the solid strip can be
            # thinner than the numerical noise in t_rib (which comes from
            # an idealized infinite centerline, not the triangulated
            # mesh), so a strict "t_rib must be inside" check spuriously
            # rejects perfectly valid slices. Distance is 0 whenever t_rib
            # genuinely falls inside a piece, so this is a strict superset
            # of the old behavior for the common case, and only kicks in
            # to rescue the near-miss / thin-strip case.
            t_min_solid = None
            t_max_solid = None
            best_dist = None
            for idxs in groups:
                piece_pts_uv = np.concatenate(
                    [pts_a_uv[idxs], pts_b_uv[idxs]], axis=0
                )
                piece_ts = np.array([t_of(uv) for uv in piece_pts_uv])
                p_t_min, p_t_max = piece_ts.min(), piece_ts.max()
                if t_rib < p_t_min:
                    dist = p_t_min - t_rib
                elif t_rib > p_t_max:
                    dist = t_rib - p_t_max
                else:
                    dist = 0.0
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    t_min_solid, t_max_solid = p_t_min, p_t_max

            if t_min_solid is None:
                continue

            # Apply margin
            t_start = t_min_solid + hole_margin
            t_end = t_max_solid - hole_margin
            if t_start >= t_end:
                continue

            t_center = (t_start + t_end) / 2.0
            available_width = t_end - t_start

            x = (d - d_min) / (d_max - d_min)
            factor = point_condition(x)

            half_w = factor * available_width / 2.0
            if half_w < 1e-6:
                continue

            hole_start = max(t_center - half_w, t_start)
            hole_end = min(t_center + half_w, t_end)
            if hole_end - hole_start < 1e-6:
                continue

            def t_to_uv(t):
                return P_uv + (t - t_rib) * bridge_dir_2d

            def uv_to_3d(uv):
                return plane_origin + uv[0] * u_ax + uv[1] * v_ax

            pt_left = uv_to_3d(t_to_uv(hole_start))
            pt_right = uv_to_3d(t_to_uv(hole_end))

            ribbons.append((d, pt_left, pt_right))

        if len(ribbons) >= 2:
            ribbons.sort(key=lambda x: x[0])
            for k in range(len(ribbons) - 1):
                _, L1, R1 = ribbons[k]
                _, L2, R2 = ribbons[k + 1]
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

    print(f"Hole mesh (analytical): {len(hole_mesh.vertices)} verts, {len(hole_mesh.faces)} faces")

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
