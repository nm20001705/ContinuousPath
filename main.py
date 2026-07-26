#!/usr/bin/env python3
import numpy as np
import trimesh
import math
import time
import os
from functools import wraps

# ============================================================
# PROFILING
# ============================================================
def timed(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        print(f"⏱️ {func.__name__} took {elapsed:.2f} seconds")
        return result
    return wrapper

# ============================================================
# PARAMETERS
# ============================================================
PLANE_DEFS = {
    'XY': (np.array([0,0,1]), np.array([1,0,0]), np.array([0,1,0])),
    'XZ': (np.array([0,1,0]), np.array([1,0,0]), np.array([0,0,1])),
    'YZ': (np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])),
}

class LWInfillParams:
    def __init__(self,
                 input_stl_path="",
                 output_stl_path="",
                 rib_spacing=20.0,
                 xy_rib_width=0.17,
                 rib_angle=30.0,
                 bridge_height=0.5,
                 bridge_width=0.5,
                 z_step=1.0,
                 cutout_margin=1.0,
                 grid_orientation=0.0,
                 primary_dir=None,
                 construction_plane='XZ',
                 rib_width=None,
                 create_holes=True):
        self.input_stl_path = input_stl_path
        self.output_stl_path = output_stl_path
        self.rib_spacing = rib_spacing
        self.rib_angle = rib_angle
        self.bridge_height = bridge_height
        self.bridge_width = bridge_width
        self.z_step = z_step
        self.cutout_margin = cutout_margin
        self.grid_orientation = grid_orientation
        self.construction_plane = construction_plane
        self.plane_normal, self.plane_axis_u, self.plane_axis_v = PLANE_DEFS[construction_plane]
        if xy_rib_width:
            self.rib_width = xy_rib_width / math.sin(math.radians(90 - rib_angle))
        else:
            self.rib_width = rib_width or 0.2
        self.create_holes = create_holes
        if primary_dir is not None:
            self.primary_dir = np.array(primary_dir)
        else:
            self.primary_dir = None
        self.show = False

# ============================================================
# ANALYTICAL GRID LINES
# ============================================================
def create_angled_grid_lines(bb, params):
    n = params.plane_normal
    au = params.plane_axis_u
    av = params.plane_axis_v

    corners = [
        np.array([bb[0,0], bb[0,1], bb[0,2]]),
        np.array([bb[1,0], bb[0,1], bb[0,2]]),
        np.array([bb[0,0], bb[1,1], bb[0,2]]),
        np.array([bb[1,0], bb[1,1], bb[0,2]]),
        np.array([bb[0,0], bb[0,1], bb[1,2]]),
        np.array([bb[1,0], bb[0,1], bb[1,2]]),
        np.array([bb[0,0], bb[1,1], bb[1,2]]),
        np.array([bb[1,0], bb[1,1], bb[1,2]]),
    ]

    def dot(v, a): return np.dot(v, a)
    u_vals = [dot(c, au) for c in corners]
    v_vals = [dot(c, av) for c in corners]
    span_u = max(u_vals) - min(u_vals)
    span_v = max(v_vals) - min(v_vals)

    pd = params.primary_dir if params.primary_dir is not None else (au if span_u >= span_v else av)
    perp_pd = np.cross(n, pd)
    perp_pd = perp_pd / np.linalg.norm(perp_pd)

    ang = math.radians(params.rib_angle)
    rot = math.radians(params.grid_orientation)

    def make_family_dir(sign):
        d = pd * math.cos(ang) + sign * perp_pd * math.sin(ang)
        k = n
        theta = rot
        d_rot = d * math.cos(theta) + np.cross(k, d) * math.sin(theta) + k * np.dot(k, d) * (1 - math.cos(theta))
        return d_rot

    d1 = make_family_dir(1)
    d2 = make_family_dir(-1)

    centre = np.array([(bb[0,0]+bb[1,0])/2, (bb[0,1]+bb[1,1])/2, (bb[0,2]+bb[1,2])/2])
    line_len = np.linalg.norm(bb[1,:] - bb[0,:]) * 2

    def generate_lines(d):
        stacking = np.cross(d, n)
        stacking = stacking / np.linalg.norm(stacking)
        proj_vals = [np.dot(c, stacking) for c in corners]
        min_p = min(proj_vals)
        max_p = max(proj_vals)
        num = int((max_p - min_p) / params.rib_spacing) + 3
        center_proj = np.dot(centre, stacking)

        lines = []
        for i in range(-1, num+1):
            offset = min_p + i * params.rib_spacing
            shift = offset - center_proj
            p0 = centre + stacking * shift
            start = p0 - d * line_len
            end   = p0 + d * line_len
            lines.append((start, end))
        return lines

    return generate_lines(d1), generate_lines(d2)

# ============================================================
# CREATE RIB MESH
# ============================================================
def create_rib_mesh(start, end, plane_normal, rib_width, y_extent):
    d = end - start
    length = np.linalg.norm(d)
    if length < 1e-8:
        return None
    d = d / length
    n_rib = np.cross(d, plane_normal)
    n_rib = n_rib / np.linalg.norm(n_rib)
    box = trimesh.creation.box(extents=[length, rib_width, y_extent])
    rot = np.eye(3)
    rot[:,0] = d
    rot[:,1] = n_rib
    rot[:,2] = plane_normal
    centre = (start + end) / 2.0
    T = np.eye(4)
    T[:3,:3] = rot
    T[:3,3] = centre - rot @ np.array([0,0,0])
    box.apply_transform(T)
    return box

# ============================================================
# EXACT FREECAD CUTOUT ALGORITHM (pure Python, trimesh)
# ============================================================
def create_rectangular_cutout_from_boundary_mesh(boundary_pts, rib_width, plane_normal, margin, tol=1e-4):
    """
    Replicates FreeCAD's create_rectangular_cutout_from_boundary.
    boundary_pts: list of 3D points (numpy arrays) on the rib plane, in order (closed loop).
    Returns a trimesh.Trimesh (the cutout face) or None.
    """
    if len(boundary_pts) < 4:
        return None

    pts = np.array(boundary_pts)

    # 1. Group by Z to find min and max
    z_vals = pts[:,2]
    z_min = np.min(z_vals)
    z_max = np.max(z_vals)

    z_tol = 1e-4
    low_mask = np.abs(z_vals - z_min) <= z_tol
    high_mask = np.abs(z_vals - z_max) <= z_tol
    low_pts = pts[low_mask]
    high_pts = pts[high_mask]

    if len(low_pts) == 0 or len(high_pts) == 0:
        return None

    low_cent = np.mean(low_pts, axis=0)
    high_cent = np.mean(high_pts, axis=0)

    z_mid = (low_cent[2] + high_cent[2]) / 2.0

    # 2. Find mid-height intersection points by interpolating along edges
    mid_pts = []
    n_pts = len(pts)
    for i in range(n_pts):
        p1 = pts[i]
        p2 = pts[(i+1) % n_pts]
        if abs(p1[2] - p2[2]) < 1e-6:
            continue
        if (p1[2] - z_mid) * (p2[2] - z_mid) < 0:
            t = (z_mid - p1[2]) / (p2[2] - p1[2])
            x = p1[0] + t * (p2[0] - p1[0])
            y = p1[1] + t * (p2[1] - p1[1])
            mid_pts.append(np.array([x, y, z_mid]))

    if len(mid_pts) < 2:
        return None

    # Remove duplicates
    uniq_mid = []
    for p in mid_pts:
        if not any(np.linalg.norm(p - q) < tol for q in uniq_mid):
            uniq_mid.append(p)

    if len(uniq_mid) < 2:
        return None

    # Choose leftmost and rightmost (by X coordinate)
    left_mid = min(uniq_mid, key=lambda p: p[0])
    right_mid = max(uniq_mid, key=lambda p: p[0])

    # 3. Four corners
    corners = np.array([low_cent, left_mid, high_cent, right_mid])

    # 4. Shrink inward by margin
    centroid = np.mean(corners, axis=0)

    def shift_towards(pt, center, margin):
        dir_vec = pt - center
        norm = np.linalg.norm(dir_vec)
        if norm < 1e-6:
            return pt
        return pt - dir_vec / norm * margin

    shifted = np.array([shift_towards(p, centroid, margin) for p in corners])

    # 5. Build a planar face (quad mesh)
    vertices = shifted
    faces = [[0, 1, 2, 3]]
    face_mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    face_mesh.fix_normals()

    # 6. Check that the face is inside the original boundary (optional)
    # We could sample points and check containment, but we trust the algorithm.

    return face_mesh

def create_prism_from_face(face_mesh, normal, height):
    """
    Extrude a planar face along the normal by height to create a solid prism.
    """
    # Use trimesh.creation.extrude_triangulation? But we have a quad face.
    # We can build the bottom and top faces, then side faces.
    vertices = face_mesh.vertices
    faces = face_mesh.faces
    # Ensure triangulated
    if len(faces[0]) == 4:
        faces = [faces[0][:3], [faces[0][0], faces[0][2], faces[0][3]]]
    # Create bottom and top
    bottom = vertices
    top = vertices + normal * height
    all_verts = np.vstack([bottom, top])
    # Faces: bottom (existing), top (reverse order), sides
    side_faces = []
    n = len(bottom)
    for i in range(n):
        j = (i+1) % n
        side_faces.append([i, j, n+j, n+i])
    # Build mesh
    mesh = trimesh.Trimesh(vertices=all_verts, faces=faces + side_faces, process=False)
    mesh.fix_normals()
    return mesh

# ============================================================
# CREATE HOLED RIB UNION (accurate, fast)
# ============================================================
@timed
def create_holed_rib_union(wing_mesh, rib_lines, params, plane_normal, rib_width, y_extent, margin, debug_dir=None):
    import time
    t0 = time.perf_counter()

    # 1. Build rib union
    t_start = time.perf_counter()
    rib_meshes = []
    for start, end in rib_lines:
        m = create_rib_mesh(start, end, plane_normal, rib_width, y_extent)
        if m is not None:
            rib_meshes.append(m)
    if not rib_meshes:
        return None
    rib_union = rib_meshes[0]
    for m in rib_meshes[1:]:
        rib_union = rib_union.union(m)
    print(f"  Rib union built in {time.perf_counter() - t_start:.2f}s")

    if not params.create_holes:
        return rib_union

    # 2. Get boundary points for each rib
    all_cutouts = []
    cutout_count = 0
    total_slice_time = 0.0
    total_cutout_face_time = 0.0
    total_prism_time = 0.0
    total_rib_loop_time = 0.0

    total_ribs = len(rib_lines)
    processed_ribs = 0
    successful_cutouts = 0

    t_loop_start = time.perf_counter()
    
    for idx, (start, end) in enumerate(rib_lines):
        t_rib_start = time.perf_counter()
        d = end - start
        length = np.linalg.norm(d)
        if length < 1e-8:
            print(f"  Rib {idx}: zero length, skipping")
            continue

        d = d / length
        n_rib = np.cross(d, plane_normal)
        n_rib = n_rib / np.linalg.norm(n_rib)
        origin = start

        # Accurate: plane-mesh intersection
        t_slice_start = time.perf_counter()
        try:
            result = trimesh.intersections.slice_mesh_plane(wing_mesh, plane_normal=n_rib, plane_origin=origin)
        except Exception as e:
            print(f"  Rib {idx}: slice_mesh_plane failed: {e}")
            continue
        t_slice = time.perf_counter() - t_slice_start
        total_slice_time += t_slice

        if result is None:
            print(f"  Rib {idx}: slice_mesh_plane returned None")
            continue

        # Extract line segments
        segments = []
        if hasattr(result, 'vertices') and hasattr(result, 'edges'):
            for edge in result.edges:
                segments.append((result.vertices[edge[0]], result.vertices[edge[1]]))
        else:
            if isinstance(result, (list, tuple)):
                for item in result:
                    if isinstance(item, (list, tuple)) and len(item) == 2:
                        segments.append((np.array(item[0]), np.array(item[1])))
            else:
                try:
                    for edge in result.edges:
                        segments.append((result.vertices[edge[0]], result.vertices[edge[1]]))
                except:
                    pass

        if not segments:
            print(f"  Rib {idx}: no boundary segments found")
            continue

        # Collect all vertices
        boundary_pts = []
        for seg in segments:
            boundary_pts.append(seg[0])
            boundary_pts.append(seg[1])

        # Remove duplicates
        unique = []
        for p in boundary_pts:
            if not any(np.linalg.norm(p - q) < 1e-3 for q in unique):
                unique.append(p)

        print(f"  Rib {idx}: segments={len(segments)}, boundary pts={len(boundary_pts)}, unique={len(unique)}")

        if len(unique) < 4:
            print(f"  Rib {idx}: not enough points for cutout (need >=4)")
            continue

        # Compute cutout face
        t_face_start = time.perf_counter()
        cutout_face = create_rectangular_cutout_from_boundary_mesh(unique, rib_width, plane_normal, margin)
        t_face = time.perf_counter() - t_face_start
        total_cutout_face_time += t_face

        if cutout_face is None:
            print(f"  Rib {idx}: cutout face creation failed")
            continue

        # Extrude to solid
        t_prism_start = time.perf_counter()
        normal = n_rib
        height = rib_width * 1.05
        prism = create_prism_from_face(cutout_face, normal, height)
        t_prism = time.perf_counter() - t_prism_start
        total_prism_time += t_prism

        if prism is not None:
            all_cutouts.append(prism)
            cutout_count += 1
            successful_cutouts += 1
            if debug_dir:
                try:
                    prism.export(os.path.join(debug_dir, f"cutout_{idx}.stl"))
                except:
                    pass
            print(f"  Rib {idx}: cutout created (slice={t_slice:.3f}s, face={t_face:.3f}s, prism={t_prism:.3f}s)")
        else:
            print(f"  Rib {idx}: prism extrusion failed")

        processed_ribs += 1
        t_rib_end = time.perf_counter()
        total_rib_loop_time += t_rib_end - t_rib_start

    print(f"\n  Summary for {processed_ribs} ribs processed:")
    print(f"    Total slice_mesh_plane time: {total_slice_time:.2f}s")
    print(f"    Total cutout face creation time: {total_cutout_face_time:.2f}s")
    print(f"    Total prism extrusion time: {total_prism_time:.2f}s")
    print(f"    Total rib loop time: {total_rib_loop_time:.2f}s")
    print(f"    Successful cutouts: {successful_cutouts} out of {len(rib_lines)} ribs")

    # 3. Fuse cutouts and subtract from rib union
    if all_cutouts:
        print("  Fusing cutouts...")
        t_fuse_start = time.perf_counter()
        cutout_union = all_cutouts[0]
        for c in all_cutouts[1:]:
            cutout_union = cutout_union.union(c)
        t_fuse = time.perf_counter() - t_fuse_start
        print(f"  Cutouts fused in {t_fuse:.2f}s")

        if debug_dir:
            try:
                cutout_union.export(os.path.join(debug_dir, "all_cutouts.stl"))
            except:
                pass

        print("  Subtracting holes from ribs...")
        t_subtract_start = time.perf_counter()
        try:
            rib_union = trimesh.boolean.difference([rib_union, cutout_union], engine='manifold')
            t_subtract = time.perf_counter() - t_subtract_start
            print(f"  Subtraction completed in {t_subtract:.2f}s")
            print("Successfully subtracted holes from ribs.")
        except Exception as e:
            print(f"Hole subtraction failed: {e}")
    else:
        print("No cutouts created.")

    print(f"  Total time for create_holed_rib_union: {time.perf_counter() - t0:.2f}s")
    return rib_union

# ============================================================
# BRIDGES (unchanged)
# ============================================================
def get_rib_boundary_z_range(wing_mesh, start, end, plane_normal, tol=1e-3):
    d = end - start
    length = np.linalg.norm(d)
    if length < 1e-6:
        return None
    d = d / length
    z_min = None
    z_max = None
    for i in range(101):
        t = i / 100.0
        pt = start + d * t * length
        if wing_mesh.contains([pt])[0]:
            if z_min is None or pt[2] < z_min:
                z_min = pt[2]
            if z_max is None or pt[2] > z_max:
                z_max = pt[2]
    if z_min is not None and z_max is not None:
        return z_min, z_max
    return None

def collect_midpoints(rib_lines, rib_z_ranges, z_step):
    data_by_rib = {}
    for idx, (start, end) in enumerate(rib_lines):
        z_min, z_max = rib_z_ranges[idx]
        if z_min is None or z_max is None:
            continue
        d = end - start
        length = np.linalg.norm(d)
        if length < 1e-8:
            continue
        d = d / length
        z_start = start[2]
        z_end = end[2]
        z_low = max(z_start, z_min) if z_start < z_end else max(z_end, z_min)
        z_high = min(z_start, z_max) if z_start > z_end else min(z_end, z_max)
        if z_high - z_low < 1e-6:
            continue
        points = []
        z = z_low
        while z <= z_high + 1e-6:
            t = (z - start[2]) / (end[2] - start[2]) if abs(end[2] - start[2]) > 1e-8 else 0
            pt = start + t * (end - start)
            points.append((z, pt))
            z += z_step
        for z_bound in (z_min, z_max):
            t = (z_bound - start[2]) / (end[2] - start[2]) if abs(end[2] - start[2]) > 1e-8 else 0
            pt = start + t * (end - start)
            if not any(abs(z_bound - p[0]) < 1e-4 for p in points):
                points.append((z_bound, pt))
        points.sort(key=lambda x: x[0])
        data_by_rib[idx] = points
    return data_by_rib

def create_bridge_meshes(data_by_rib, plane_normal, rib_width, bridge_height, y_extent):
    bridges = []
    for idx, pts in data_by_rib.items():
        if len(pts) < 2:
            continue
        for i in range(len(pts)-1):
            p1 = pts[i][1]
            p2 = pts[i+1][1]
            d = p2 - p1
            length = np.linalg.norm(d)
            if length < 1e-8:
                continue
            d = d / length
            n_rib = np.cross(d, plane_normal)
            n_rib = n_rib / np.linalg.norm(n_rib)
            width = rib_width
            box = trimesh.creation.box(extents=[length, width, bridge_height])
            rot = np.eye(3)
            rot[:,0] = d
            rot[:,1] = n_rib
            rot[:,2] = plane_normal
            centre = (p1 + p2) / 2.0
            T = np.eye(4)
            T[:3,:3] = rot
            T[:3,3] = centre - rot @ np.array([0,0,0])
            box.apply_transform(T)
            bridges.append(box)
    return bridges

# ============================================================
# MAIN
# ============================================================
@timed
def main(params):
    if not params.input_stl_path:
        raise ValueError("input_stl_path must be provided")
    wing_mesh = trimesh.load(params.input_stl_path)
    if wing_mesh is None or not isinstance(wing_mesh, trimesh.Trimesh):
        raise ValueError("Failed to load STL or loaded object is not a mesh")

    # Repair mesh (soft)
    original_mesh = wing_mesh
    try:
        if not wing_mesh.is_watertight:
            repaired = trimesh.repair.fill_holes(wing_mesh)
            if repaired is not None:
                wing_mesh = repaired
        wing_mesh = trimesh.repair.fix_normals(wing_mesh)
        wing_mesh = wing_mesh.merge_vertices()
    except Exception as e:
        print(f"Warning: mesh repair failed ({e}), using original mesh")
        wing_mesh = original_mesh

    debug_dir = os.path.dirname(params.output_stl_path)
    if debug_dir:
        debug_dir = os.path.join(debug_dir, "debug")
        os.makedirs(debug_dir, exist_ok=True)

    bb = wing_mesh.bounds
    y_min = bb[0,1]
    y_max = bb[1,1]
    y_extent = (y_max - y_min) * 1.2

    # 1. Grid lines
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_lines = lines1 + lines2

    # 2. Rib Z ranges (for bridges)
    rib_z_ranges = []
    for start, end in all_lines:
        z_range = get_rib_boundary_z_range(wing_mesh, start, end, params.plane_normal)
        if z_range is None:
            rib_z_ranges.append((None, None))
        else:
            rib_z_ranges.append(z_range)

    # 3. Holed rib union (accurate)
    rib_union = create_holed_rib_union(wing_mesh, all_lines, params, params.plane_normal,
                                       params.rib_width, y_extent, params.cutout_margin,
                                       debug_dir=debug_dir)
    if rib_union is None:
        print("No ribs generated – exiting.")
        return

    # 4. Subtract ribs from wing
    print("Subtracting ribs from wing...")
    try:
        cut_wing = trimesh.boolean.difference([wing_mesh, rib_union], engine='manifold')
    except Exception as e:
        print(f"Boolean difference failed: {e}")
        cut_wing = wing_mesh

    # 5. Bridges
    data_by_rib = collect_midpoints(all_lines, rib_z_ranges, params.z_step)
    bridge_meshes = create_bridge_meshes(data_by_rib, params.plane_normal,
                                         params.rib_width, params.bridge_height, y_extent)

    if bridge_meshes:
        print("Adding bridges...")
        bridge_union = bridge_meshes[0]
        for m in bridge_meshes[1:]:
            bridge_union = bridge_union.union(m)
        try:
            final_mesh = trimesh.boolean.union([cut_wing, bridge_union], engine='manifold')
        except Exception as e:
            print(f"Union failed: {e}")
            final_mesh = cut_wing
    else:
        final_mesh = cut_wing

    # 6. Export
    if params.output_stl_path:
        final_mesh.export(params.output_stl_path)
        print(f"Exported STL to {params.output_stl_path}")

    if params.show:
        final_mesh.show()

# ============================================================
# RUN
# ============================================================
if __name__ == "__main__":
    params = LWInfillParams(
        input_stl_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\files\Fin.stl",
        output_stl_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin_out.stl",
        rib_spacing=20.0,
        xy_rib_width=0.17,
        rib_angle=30.0,
        bridge_height=0.5,
        bridge_width=0.5,
        z_step=1.0,
        cutout_margin=1.0,
        construction_plane='XZ',
        create_holes=True,
    )
    main(params)