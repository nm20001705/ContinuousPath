#!/usr/bin/env python3
import numpy as np
import trimesh
import math
import time
import os
import tempfile
from functools import wraps

# ============================================================
# STEP EXPORTER (using cadquery)
# ============================================================
def export_to_step(mesh, path):
    """
    Export a trimesh mesh as a STEP file using cadquery.
    Saves mesh to a temporary STL, imports it with cadquery, and exports to STEP.
    Falls back to STL if cadquery is not available or import fails.
    """
    try:
        import cadquery as cq
        from cadquery import importers
        # Write mesh to a temporary STL file
        with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp:
            tmp_stl = tmp.name
        mesh.export(tmp_stl)
        # Import STL using cadquery's importers
        try:
            # Newer versions: importers.importMesh
            cq_object = importers.importMesh(tmp_stl)
        except AttributeError:
            # Older versions: importers.importStl
            try:
                cq_object = importers.importStl(tmp_stl)
            except AttributeError:
                # If neither works, fallback
                raise ImportError("No suitable import function found in cadquery")
        # Export to STEP
        cq.exporters.export(cq_object, path)
        os.unlink(tmp_stl)
        print(f"  Exported STEP to {path}")
    except Exception as e:
        print(f"  STEP export failed: {e}")
        # Fallback to STL with .stl extension
        stl_path = path.rsplit('.', 1)[0] + '.stl'
        mesh.export(stl_path)
        print(f"  Saved STL as {stl_path}")

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
# CUTOUT FUNCTIONS
# ============================================================
def create_rectangular_cutout_from_boundary_mesh(boundary_pts, rib_width, plane_normal, margin, tol=1e-4):
    if len(boundary_pts) < 4:
        return None
    pts = np.array(boundary_pts)
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
    uniq_mid = []
    for p in mid_pts:
        if not any(np.linalg.norm(p - q) < tol for q in uniq_mid):
            uniq_mid.append(p)
    if len(uniq_mid) < 2:
        return None
    left_mid = min(uniq_mid, key=lambda p: p[0])
    right_mid = max(uniq_mid, key=lambda p: p[0])
    corners = np.array([low_cent, left_mid, high_cent, right_mid])
    centroid = np.mean(corners, axis=0)
    def shift_towards(pt, center, margin):
        dir_vec = pt - center
        norm = np.linalg.norm(dir_vec)
        if norm < 1e-6:
            return pt
        return pt - dir_vec / norm * margin
    shifted = np.array([shift_towards(p, centroid, margin) for p in corners])
    vertices = shifted
    faces = [[0, 1, 2, 3]]
    face_mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    face_mesh.fix_normals()
    return face_mesh

def create_prism_from_face(face_mesh, normal, height):
    """
    Extrude a planar face along its normal to create a watertight prism.
    """
    # Ensure face is planar (we assume it is)
    vertices = face_mesh.vertices
    if len(vertices) < 3:
        return None

    # Compute a local 2D coordinate system
    origin = vertices[0]
    # Use first edge as u‑axis
    u = vertices[1] - origin
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    if np.linalg.norm(v) < 1e-8:
        v = np.cross(normal, [1,0,0])
        if np.linalg.norm(v) < 1e-8:
            v = np.cross(normal, [0,1,0])
        v = v / np.linalg.norm(v)
        u = np.cross(v, normal)
        u = u / np.linalg.norm(u)

    # Project vertices to 2D
    coords_2d = []
    for p in vertices:
        dx = np.dot(p - origin, u)
        dy = np.dot(p - origin, v)
        coords_2d.append([dx, dy])
    coords_2d = np.array(coords_2d)

    # Order vertices to form a convex polygon (or use Shapely to handle ordering)
    from shapely.geometry import Polygon as ShapelyPolygon
    poly = ShapelyPolygon(coords_2d)
    if not poly.is_valid:
        poly = poly.buffer(0)
    ext_vertices = np.array(poly.exterior.coords)[:-1]  # remove repeated last point

    # Map back to 3D
    polygon_3d = []
    for pt2d in ext_vertices:
        x, y = pt2d
        p3d = origin + x * u + y * v
        polygon_3d.append(p3d)

    # Create a path and extrude
    path = trimesh.path.Path3D(
        entities=[trimesh.path.entities.Line(np.arange(len(polygon_3d)+1))],
        vertices=np.vstack([polygon_3d, polygon_3d[0]])
    )
    try:
        extruded = trimesh.creation.extrude_polygon(path, height, transform=None)
    except:
        # Fallback: convex hull of the extruded vertices
        bottom = polygon_3d
        top = [p + normal * height for p in bottom]
        all_verts = np.vstack([bottom, top])
        extruded = trimesh.convex.convex_hull(all_verts)
    return extruded

# ============================================================
# CREATE HOLED RIB UNION (with STEP debug export)
# ============================================================
@timed
def create_holed_rib_union(wing_mesh, rib_lines, params, plane_normal, rib_width, y_extent, margin, debug_dir=None, debug_format='step'):
    t0 = time.perf_counter()

    # ---- Simplify wing mesh for slicing ----
    print("  Simplifying wing mesh for slicing...")
    t_simp = time.perf_counter()
    wing_mesh_simple = wing_mesh.simplify_quadric_decimation(face_count=20000)
    print(f"  Simplified mesh from {len(wing_mesh.faces)} to {len(wing_mesh_simple.faces)} faces in {time.perf_counter()-t_simp:.2f}s")

    # ---- STEP 1: Build all rib meshes ----
    print("  Building rib meshes...")
    t_start = time.perf_counter()
    rib_meshes = []
    for idx, (start, end) in enumerate(rib_lines):
        m = create_rib_mesh(start, end, plane_normal, rib_width, y_extent)
        if m is not None:
            rib_meshes.append(m)
        if (idx + 1) % 10 == 0:
            print(f"    Built {idx + 1} rib meshes")
    if not rib_meshes:
        return None
    print(f"  Built {len(rib_meshes)} rib meshes in {time.perf_counter() - t_start:.2f}s")

    # ---- STEP 2: Union all rib meshes ----
    print("  Unioning rib meshes...")
    t_union = time.perf_counter()
    rib_union = rib_meshes[0]
    for m in rib_meshes[1:]:
        rib_union = rib_union.union(m)
    print(f"  Rib union built in {time.perf_counter() - t_union:.2f}s")

    if not params.create_holes:
        return rib_union

    # ---- STEP 3: Process each rib to create cutouts and subtract individually ----
    cutout_count = 0
    print(f"  Processing cutouts, number of rib_lines = {len(rib_lines)}")
    t_cutout_start = time.perf_counter()

    for rib_idx, (start, end) in enumerate(rib_lines):
        d = end - start
        length = np.linalg.norm(d)
        if length < 1e-8:
            continue
        d = d / length
        n_rib = np.cross(d, plane_normal)
        n_rib = n_rib / np.linalg.norm(n_rib)
        origin = start

        # 3a: Plane-mesh intersection
        t_slice = time.perf_counter()
        try:
            result = trimesh.intersections.slice_mesh_plane(wing_mesh_simple, plane_normal=n_rib, plane_origin=origin)
        except Exception as e:
            print(f"  Rib {rib_idx}: slice_mesh_plane failed: {e}")
            continue
        t_slice = time.perf_counter() - t_slice
        if t_slice > 0.1:
            print(f"  Rib {rib_idx}: slice_mesh_plane took {t_slice:.2f}s")

        if result is None:
            continue

        # Extract segments
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
            print(f"  Rib {rib_idx}: no boundary segments found")
            continue

        # Collect vertices
        boundary_pts = []
        for seg in segments:
            boundary_pts.append(seg[0])
            boundary_pts.append(seg[1])

        # Remove duplicates
        pts_array = np.array(boundary_pts)
        scaled = np.round(pts_array / 1e-3).astype(np.int64)
        _, unique_indices = np.unique(scaled, axis=0, return_index=True)
        unique = pts_array[np.sort(unique_indices)]

        if len(unique) < 4:
            print(f"  Rib {rib_idx}: not enough points for cutout (need >=4)")
            continue

        # 3b: Create cutout face
        t_face = time.perf_counter()
        cutout_face = create_rectangular_cutout_from_boundary_mesh(unique, rib_width, plane_normal, margin)
        t_face = time.perf_counter() - t_face
        if t_face > 0.1:
            print(f"  Rib {rib_idx}: face creation took {t_face:.2f}s")
        if cutout_face is None:
            print(f"  Rib {rib_idx}: cutout face creation failed")
            continue

        # 3c: Extrude to prism
        t_prism = time.perf_counter()
        normal = n_rib
        height = rib_width * 1.05
        prism = create_prism_from_face(cutout_face, normal, height)
        t_prism = time.perf_counter() - t_prism
        if t_prism > 0.1:
            print(f"  Rib {rib_idx}: prism extrusion took {t_prism:.2f}s")

        if prism is None:
            print(f"  Rib {rib_idx}: prism extrusion failed")
            continue

        # Validate volume and repair if needed
        if not prism.is_volume:
            repaired = trimesh.repair.fill_holes(prism)
            if repaired is not None and repaired.is_volume:
                prism = repaired
            else:
                # Try convex hull as last resort
                try:
                    hull = trimesh.convex.convex_hull(prism.vertices)
                    if hull.is_volume:
                        prism = hull
                    else:
                        print(f"  Rib {rib_idx}: cutout prism cannot be made a volume – skipping")
                        continue
                except:
                    print(f"  Rib {rib_idx}: cutout prism cannot be made a volume – skipping")
                    continue

        # Now we have a volume, subtract it from rib_union directly
        print(f"  Rib {rib_idx}: cutout created, subtracting from rib union...")
        try:
            # Try manifold engine first
            rib_union = trimesh.boolean.difference([rib_union, prism], engine='manifold')
            cutout_count += 1
            if debug_dir:
                ext = '.step' if debug_format.lower() == 'step' else '.stl'
                try:
                    if debug_format.lower() == 'step':
                        export_to_step(prism, os.path.join(debug_dir, f"cutout_{rib_idx}.step"))
                    else:
                        prism.export(os.path.join(debug_dir, f"cutout_{rib_idx}.stl"))
                except Exception as e:
                    print(f"    Failed to export cutout {rib_idx}: {e}")
        except Exception as e:
            print(f"  Rib {rib_idx}: subtraction failed: {e}")

    print(f"  Cutout processing took {time.perf_counter() - t_cutout_start:.2f}s")
    print(f"Total cutouts subtracted: {cutout_count}")

    # ---- Final repair of rib_union to ensure it's a volume ----
    if not rib_union.is_volume:
        print("  Repairing rib union to make it a volume...")
        repaired = trimesh.repair.fill_holes(rib_union)
        if repaired is not None and repaired.is_volume:
            rib_union = repaired
        else:
            # Try convex hull (may change shape drastically)
            try:
                hull = trimesh.convex.convex_hull(rib_union.vertices)
                if hull.is_volume:
                    rib_union = hull
                else:
                    print("  Warning: rib union could not be made a volume.")
            except:
                print("  Warning: rib union could not be made a volume.")

    if debug_dir:
        # Export final rib_union with holes
        try:
            if debug_format.lower() == 'step':
                export_to_step(rib_union, os.path.join(debug_dir, "rib_union_with_holes.step"))
            else:
                rib_union.export(os.path.join(debug_dir, "rib_union_with_holes.stl"))
        except Exception as e:
            print(f"  Failed to export final rib union: {e}")

    print(f"  Total time for create_holed_rib_union: {time.perf_counter() - t0:.2f}s")
    return rib_union

# ============================================================
# BRIDGES
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

    # 3. Holed rib union (with STEP debug export)
    rib_union = create_holed_rib_union(wing_mesh, all_lines, params, params.plane_normal,
                                       params.rib_width, y_extent, params.cutout_margin,
                                       debug_dir=debug_dir, debug_format='step')
    if rib_union is None:
        print("No ribs generated – exiting.")
        return

    # 4. Subtract ribs from wing (only if rib_union is a volume and wing_mesh is watertight)
    print("Subtracting ribs from wing...")
    if rib_union.is_volume and wing_mesh.is_watertight:
        try:
            cut_wing = trimesh.boolean.difference([wing_mesh, rib_union], engine='manifold')
        except Exception as e:
            print(f"Boolean difference failed: {e}")
            cut_wing = wing_mesh
    else:
        print("Skipping boolean difference: rib union or wing mesh is not watertight.")
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
        # Ensure bridge union is a volume
        if not bridge_union.is_volume:
            repaired = trimesh.repair.fill_holes(bridge_union)
            if repaired is not None and repaired.is_volume:
                bridge_union = repaired
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
        input_stl_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\WingR1_msv_orient.stl",
        output_stl_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\WingR1_msv_out_orient.stl",
        rib_spacing=20.0,
        xy_rib_width=0.17,
        rib_angle=30.0,
        bridge_height=0.5,
        bridge_width=0.5,
        z_step=1.0,
        cutout_margin=0.0,
        construction_plane='XZ',
        create_holes=True,
    )
    main(params)