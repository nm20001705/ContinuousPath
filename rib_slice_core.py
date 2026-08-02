# rib_slice_core.py
#
# Shared low-level machinery for slicing individual rib-segment meshes at
# each Z step, clustering the raw intersection edges into connected
# pieces, and picking the piece nearest the rib centerline (t_rib).
#
# Used by BOTH hole_utils_analytical.py and bridge_utils_analytical.py so
# the two stay in sync on how voids are handled -- they only differ in
# how wide a strip they cut out of the chosen piece.

import numpy as np
import trimesh


def basis_vectors(primary_dir_np):
    """Global (u,v) basis perpendicular to prim. Must match the basis
    used elsewhere in the pipeline (precompute_slices etc.)."""
    prim = primary_dir_np / np.linalg.norm(primary_dir_np)
    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)
    return prim, u_ax, v_ax


def cluster_segments(pts_a, pts_b, tol=5):
    """
    Group segment indices into connected components based on shared
    endpoints. pts_a/pts_b are (N,2) arrays of uv endpoints for N raw
    mesh_plane segments. Returns a list of index-arrays, one per piece.
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


def prepare_segment(seg, prim, u_ax, v_ax):
    """Precompute p0, dir_rib, bridge_dir_2d (2D) AND the raw 3D frame
    (line_dir_3d, slab_normal, plane_offset) for one rib segment dict.
    Returns None if degenerate. The 3D frame is what solidify_utils
    needs to extrude this rib line's strip into a real volume."""
    p0 = seg['p0']
    dir_rib = seg['dir']
    slab_normal = seg['slab_normal']

    line_dir_3d = np.cross(slab_normal, prim)
    norm_l = np.linalg.norm(line_dir_3d)
    if norm_l < 1e-8:
        return None
    line_dir_3d /= norm_l
    bridge_dir_2d = np.array([np.dot(line_dir_3d, u_ax), np.dot(line_dir_3d, v_ax)])
    bd_norm = np.linalg.norm(bridge_dir_2d)
    if bd_norm < 1e-8:
        return None
    bridge_dir_2d /= bd_norm

    plane_offset = np.dot(p0, slab_normal)

    return {
        'p0': p0,
        'dir_rib': dir_rib,
        'bridge_dir_2d': bridge_dir_2d,
        'line_dir_3d': line_dir_3d,
        'slab_normal': slab_normal,
        'plane_offset': plane_offset,
    }


def uv_of(pt3d, d, prim, u_ax, v_ax):
    delta = pt3d - d * prim
    return np.array([np.dot(delta, u_ax), np.dot(delta, v_ax)])


def t_to_3d(t, P_uv, t_rib, bridge_dir_2d, plane_origin, u_ax, v_ax):
    """Reconstruct a 3D point at parameter t along the straight bridge
    line through the rib centerline's projection."""
    uv = P_uv + (t - t_rib) * bridge_dir_2d
    return plane_origin + uv[0] * u_ax + uv[1] * v_ax


def iter_solid_pieces(seg_mesh, seg, prim, u_ax, v_ax, z_vals):
    """
    Generator over Z-slices of a single rib-segment mesh. For each slice
    where mesh_plane returns geometry, clusters the raw edges into
    connected pieces and picks the piece CLOSEST to the rib centerline's
    projection (t_rib) -- distance 0 if t_rib genuinely falls inside it.
    Nearest-piece (rather than strict containment) matters near thin
    sections like the trailing edge, where t_rib (from an idealized
    infinite centerline) can land a hair outside the actual triangulated
    piece due to numerical noise.

    A slice is skipped only if mesh_plane returns nothing to cluster at
    all (e.g. this rib segment simply doesn't span that far in Z).

    Yields dicts with keys:
        d, plane_origin, t_min_solid, t_max_solid,
        t_rib, P_uv, bridge_dir_2d, d_min, d_max
    """
    if seg_mesh is None or len(seg_mesh.vertices) == 0:
        return

    prep = prepare_segment(seg, prim, u_ax, v_ax)
    if prep is None:
        return
    p0 = prep['p0']
    dir_rib = prep['dir_rib']
    bridge_dir_2d = prep['bridge_dir_2d']
    line_dir_3d = prep['line_dir_3d']
    slab_normal = prep['slab_normal']
    plane_offset = prep['plane_offset']

    def t_of(uv):
        return np.dot(uv, bridge_dir_2d)

    proj = seg_mesh.vertices @ prim
    d_min = proj.min()
    d_max = proj.max()
    if d_max - d_min < 1e-8:
        return

    start_idx = np.searchsorted(z_vals, d_min)
    end_idx = np.searchsorted(z_vals, d_max, side='right')

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

        pts_a_uv = np.array([uv_of(p, d, prim, u_ax, v_ax) for p in lines_arr[:, 0, :]])
        pts_b_uv = np.array([uv_of(p, d, prim, u_ax, v_ax) for p in lines_arr[:, 1, :]])

        denom = np.dot(prim, dir_rib)
        if abs(denom) < 1e-8:
            continue
        t_plane = np.dot(prim, plane_origin - p0) / denom
        P_3d = p0 + t_plane * dir_rib
        P_uv = uv_of(P_3d, d, prim, u_ax, v_ax)
        t_rib = t_of(P_uv)

        axis_d = dir_rib / denom
        origin_const = p0 - (np.dot(prim, p0) / denom) * dir_rib

        groups = cluster_segments(pts_a_uv, pts_b_uv)

        t_min_solid = None
        t_max_solid = None
        best_dist = None
        for idxs in groups:
            piece_pts_uv = np.concatenate([pts_a_uv[idxs], pts_b_uv[idxs]], axis=0)
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

        yield {
            'd': d,
            'plane_origin': plane_origin,
            't_min_solid': t_min_solid,
            't_max_solid': t_max_solid,
            't_rib': t_rib,
            'P_uv': P_uv,
            'bridge_dir_2d': bridge_dir_2d,
            'd_min': d_min,
            'd_max': d_max,
            'line_dir_3d': line_dir_3d,
            'slab_normal': slab_normal,
            'plane_offset': plane_offset,
            'axis_d': axis_d,
            'origin_const': origin_const
        }

def hole_width_interval(piece, point_condition, hole_margin, thickness):
    d = piece['d']
    d_min = piece['d_min']
    d_max = piece['d_max']

    eff_margin = hole_margin + 0.5 * thickness

    # d-direction margin: real clearance from THIS rib's own crossing
    # surface at top/bottom, instead of relying on point_condition to
    # taper to zero exactly at the raw geometric boundary.
    d_start = d_min + eff_margin
    d_end = d_max - eff_margin
    if d_start >= d_end or d < d_start or d > d_end:
        return None

    # t-direction margin: same eff_margin, from the neighboring
    # crossing rib's solid surface (unchanged from before).
    t_start = piece['t_min_solid'] + eff_margin
    t_end = piece['t_max_solid'] - eff_margin
    if t_start >= t_end:
        return None

    t_center = (t_start + t_end) / 2.0
    available_width = t_end - t_start

    x = (d - d_start) / (d_end - d_start) if d_end - d_start > 1e-8 else 0.5
    factor = point_condition(x)

    half_w = factor * available_width / 2.0
    if half_w < 1e-6:
        return None

    hs = max(t_center - half_w, t_start)
    he = min(t_center + half_w, t_end)
    if he - hs < 1e-6:
        return None
    return hs, he

def bridge_width_interval(piece, bridge_height):
    """Width policy for bridges: constant bridge_height, centered in the
    piece, using the piece's raw solid extent. No margin term -- unlike
    holes, bridges are the load-carrying material itself, so they get no
    clearance inset. Returns None if it doesn't fit (no clamping)."""
    t_start = piece['t_min_solid']
    t_end = piece['t_max_solid']
    if t_end - t_start < bridge_height:
        return None
    t_center = (t_start + t_end) / 2.0
    return t_center - bridge_height / 2.0, t_center + bridge_height / 2.0


def collect_line_intervals(seg_mesh, seg, prim, u_ax, v_ax, z_vals, width_policy_fn):
    """
    Runs iter_solid_pieces for one rib segment and applies width_policy_fn
    to each slice. width_policy_fn takes ONLY `piece` -- callers close
    over thickness/margin/etc. themselves (see hole_utils.py /
    bridge_utils.py). Returns (intervals, frame) where:
      intervals : list of (d0, d1, t_start, t_end) -- consecutive Z-steps
                  paired up, ready for solidify_rib_line. Note: this uses
                  the axis-aligned bounding rectangle of each pair's four
                  t-values rather than the exact trapezoid, which is a
                  deliberate (slightly conservative) simplification --
                  shapely's box() is axis-aligned, and the small amount
                  of extra area at a taper is negligible next to
                  hole_margin/bridge_height scale features.
      frame     : dict with plane_offset, prim, line_dir_3d, slab_normal
                  (None if no valid slices were found).
    """
    rows = []  # (d, t_start, t_end)
    frame = None

    for piece in iter_solid_pieces(seg_mesh, seg, prim, u_ax, v_ax, z_vals):
            result = width_policy_fn(piece)
            if result is None:
                continue
            t_start, t_end = result
            t_rib = piece['t_rib']
            # Store relative to this slice's centerline-crossing point.
            # origin_const + d*axis_d already reconstructs P_3d(d) exactly
            # (the point where the rib centerline crosses height d); the
            # remaining in-wall offset is (t - t_rib(d)) along line_dir_3d,
            # NOT t itself. t_rib(d) is generally nonzero and drifts with d,
            # so using raw t here was adding a d-dependent shift to every
            # strip -- this is what was showing up as bridges/holes drifting
            # away from their correct position as height increased.
            rows.append((piece['d'], t_start - t_rib, t_end - t_rib))
            if frame is None:
                frame = {
                    'origin_const': piece['origin_const'],
                    'axis_d': piece['axis_d'],
                    'line_dir_3d': piece['line_dir_3d'],
                    'slab_normal': piece['slab_normal'],
                }

    if len(rows) < 2:
        return [], frame

    rows.sort(key=lambda r: r[0])
    intervals = []
    for k in range(len(rows) - 1):
        d0, t0a, t0b = rows[k]
        d1, t1a, t1b = rows[k + 1]
        t_lo = min(t0a, t0b, t1a, t1b)
        t_hi = max(t0a, t0b, t1a, t1b)
        intervals.append((d0, d1, t_lo, t_hi))

    return intervals, frame


def ribbons_to_mesh(ribbons):
    """
    ribbons: list of (d, pt_left, pt_right) tuples, sorted by d, for a
    SINGLE rib segment. Returns (vertices, faces) with 0-based LOCAL face
    indices -- caller offsets them into the global vertex buffer.
    """
    verts = []
    faces = []
    offset = 0
    for k in range(len(ribbons) - 1):
        _, L1, R1 = ribbons[k]
        _, L2, R2 = ribbons[k + 1]
        v0, v1, v2, v3 = offset, offset + 1, offset + 2, offset + 3
        verts.extend([L1, R1, R2, L2])
        faces.append([v0, v1, v2])
        faces.append([v0, v2, v3])
        offset += 4
    return verts, faces