#!/usr/bin/env python3
import numpy as np
import trimesh
import math
import time
import os
import tempfile
from functools import wraps
import shapely.geometry as sg
from shapely.ops import polygonize

# ============================================================
# MESH REPAIR
# ============================================================
def repair_mesh(mesh):
    if mesh is None:
        return None
    if isinstance(mesh, trimesh.Scene):
        meshes = [g for g in mesh.geometry.values() if isinstance(g, trimesh.Trimesh)]
        if not meshes:
            return None
        if len(meshes) == 1:
            mesh = meshes[0]
        else:
            mesh = trimesh.util.concatenate(meshes)
            if mesh is None:
                return None
    if not isinstance(mesh, trimesh.Trimesh):
        return None
    if len(mesh.faces) == 0:
        return None
    try:
        merged = mesh.merge_vertices()
        if merged is not None:
            mesh = merged
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
    for max_hole in [0.01, 0.05, 0.1, 0.5, 1.0]:
        try:
            repaired = trimesh.repair.fill_holes(mesh, max_hole=max_hole)
            if repaired is not None and repaired.is_watertight:
                mesh = repaired
                break
        except:
            continue
    try:
        trimesh.repair.wind_watertight(mesh)
    except:
        pass
    try:
        trimesh.repair.broken_faces(mesh)
    except:
        pass
    try:
        mesh.fix_normals()
    except:
        pass
    if not mesh.is_watertight:
        try:
            hull = trimesh.convex.convex_hull(mesh.vertices)
            if hull.is_watertight:
                mesh = hull
        except:
            pass
    return mesh

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
# CREATE HOLED RIB UNION (with robust cutouts using Shapely)
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
        try:
            rib_union = trimesh.boolean.union([rib_union, m], engine='manifold')
        except:
            try:
                rib_union = trimesh.boolean.union([rib_union, m], engine='scad')
            except:
                print("  Skipping a rib due to union failure")
                continue
    print(f"  Rib union built in {time.perf_counter() - t_union:.2f}s")

    if not params.create_holes:
        return rib_union

    # ---- STEP 3: Collect cutout prisms for a single subtraction ----
    cutout_prisms = []
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

        # Plane-mesh intersection
        try:
            path = trimesh.intersections.slice_mesh_plane(wing_mesh_simple, plane_normal=n_rib, plane_origin=origin)
        except Exception as e:
            print(f"  Rib {rib_idx}: slice_mesh_plane failed: {e}")
            continue

        if path is None or len(path.vertices) < 3:
            continue

        # Build local frame for the rib plane
        # Use the first edge as u-axis
        verts = path.vertices
        edges = path.edges
        # Find a suitable u vector from the first edge
        u = verts[edges[0][1]] - verts[edges[0][0]]
        u_norm = np.linalg.norm(u)
        if u_norm < 1e-8:
            # fallback
            u = np.cross(n_rib, [1,0,0])
            if np.linalg.norm(u) < 1e-8:
                u = np.cross(n_rib, [0,1,0])
            u = u / np.linalg.norm(u)
        else:
            u = u / u_norm
        v = np.cross(n_rib, u)
        if np.linalg.norm(v) < 1e-8:
            v = np.cross(n_rib, [0,1,0])
            v = v / np.linalg.norm(v)
            u = np.cross(v, n_rib)
            u = u / np.linalg.norm(u)
        v = v / np.linalg.norm(v)

        # Project vertices to 2D
        coords_2d = []
        for p in verts:
            dx = np.dot(p - origin, u)
            dy = np.dot(p - origin, v)
            coords_2d.append((dx, dy))

        # Build line segments for polygonize
        lines = []
        for edge in edges:
            p1 = coords_2d[edge[0]]
            p2 = coords_2d[edge[1]]
            lines.append(sg.LineString([p1, p2]))
        merged = sg.MultiLineString(lines)
        polygons = list(polygonize(merged))
        if not polygons:
            continue
        # Take the largest polygon
        poly = max(polygons, key=lambda p: p.area)
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty or poly.area < 1e-8:
            continue

        # Extrude the 2D polygon along the rib normal
        height = rib_width * 1.05
        # Build transform: local Z -> normal, local origin -> origin
        transform = np.eye(4)
        transform[:3,:3] = np.column_stack([u, v, n_rib])
        transform[:3,3] = origin

        try:
            extruded = trimesh.creation.extrude_polygon(poly, height, transform=transform)
        except Exception as e:
            print(f"  Rib {rib_idx}: extrusion failed: {e}")
            continue

        if extruded is None:
            continue
        # Ensure the extruded mesh is a volume
        if not extruded.is_volume:
            # Try to repair
            repaired = trimesh.repair.fill_holes(extruded)
            if repaired is not None and repaired.is_volume:
                extruded = repaired
            else:
                # Fallback: convex hull of the union of bottom and top faces
                try:
                    # Get bottom and top vertices from the extrusion
                    # The extrusion creates a mesh, we can extract the vertices that lie on the top and bottom planes
                    # But simpler: create a convex hull of all vertices
                    hull = trimesh.convex.convex_hull(extruded.vertices)
                    if hull.is_volume:
                        extruded = hull
                    else:
                        continue
                except:
                    continue

        # Ensure normals are consistent
        extruded.fix_normals()
        if extruded.is_volume:
            cutout_prisms.append(extruded)
            if debug_dir:
                try:
                    if debug_format.lower() == 'step':
                        export_to_step(extruded, os.path.join(debug_dir, f"cutout_{rib_idx}.step"))
                    else:
                        extruded.export(os.path.join(debug_dir, f"cutout_{rib_idx}.stl"))
                except:
                    pass
        else:
            print(f"  Rib {rib_idx}: extruded mesh is not a volume")

    print(f"  Cutout processing took {time.perf_counter() - t_cutout_start:.2f}s")
    print(f"  Collected {len(cutout_prisms)} valid cutout prisms.")

    # ---- STEP 4: Union cutouts and subtract from rib_union ----
    if cutout_prisms:
        print("  Unioning cutout prisms...")
        cutout_union = cutout_prisms[0]
        for cp in cutout_prisms[1:]:
            try:
                cutout_union = trimesh.boolean.union([cutout_union, cp], engine='manifold')
            except:
                try:
                    cutout_union = trimesh.boolean.union([cutout_union, cp], engine='scad')
                except:
                    print("  Skipping a cutout union due to failure")
                    continue
        if cutout_union.is_volume and rib_union.is_volume:
            print("  Subtracting cutout union from rib union...")
            try:
                rib_union = trimesh.boolean.difference([rib_union, cutout_union], engine='manifold')
            except:
                try:
                    rib_union = trimesh.boolean.difference([rib_union, cutout_union], engine='scad')
                except Exception as e:
                    print(f"  Subtraction failed: {e}")
            # Repair after large subtraction
            if not rib_union.is_volume:
                print("  Repairing rib union after cutout subtraction...")
                rib_union = trimesh.repair.fill_holes(rib_union)
                if rib_union is None or not rib_union.is_volume:
                    try:
                        hull = trimesh.convex.convex_hull(rib_union.vertices)
                        if hull.is_volume:
                            rib_union = hull
                    except:
                        pass
        else:
            print("  Skipping cutout subtraction because cutout union or rib union is not watertight.")

    # ---- Final repair of rib_union ----
    if not rib_union.is_volume:
        print("  Final repair of rib union to make it a volume...")
        repaired = trimesh.repair.fill_holes(rib_union)
        if repaired is not None and repaired.is_volume:
            rib_union = repaired
        else:
            try:
                hull = trimesh.convex.convex_hull(rib_union.vertices)
                if hull.is_volume:
                    rib_union = hull
                else:
                    print("  Warning: rib union could not be made a volume.")
            except:
                print("  Warning: rib union could not be made a volume.")

    if debug_dir:
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

    print("Loading mesh...")
    loaded = trimesh.load(params.input_stl_path)

    # If the loaded object is a Scene, extract the first mesh
    if isinstance(loaded, trimesh.Scene):
        meshes = [g for g in loaded.geometry.values() if isinstance(g, trimesh.Trimesh)]
        if not meshes:
            raise ValueError("No Trimesh found in the Scene")
        wing_mesh = meshes[0]
    else:
        wing_mesh = loaded

    if wing_mesh is None or not isinstance(wing_mesh, trimesh.Trimesh):
        raise ValueError("Failed to load a valid Trimesh object")

    print(f"Mesh vertices: {len(wing_mesh.vertices)}, faces: {len(wing_mesh.faces)}")
    print(f"Is watertight: {wing_mesh.is_watertight}")
    print(f"Is volume: {wing_mesh.is_volume}")

    # Repair the wing mesh to make it watertight if possible
    print("Repairing wing mesh...")
    wing_mesh = repair_mesh(wing_mesh)
    if wing_mesh is None:
        raise ValueError("Mesh repair failed")
    print(f"After repair: watertight={wing_mesh.is_watertight}, volume={wing_mesh.is_volume}")

    debug_dir = os.path.dirname(params.output_stl_path)
    if debug_dir:
        debug_dir = os.path.join(debug_dir, "debug")
        os.makedirs(debug_dir, exist_ok=True)

    bb = wing_mesh.bounds
    y_min = bb[0,1]
    y_max = bb[1,1]
    y_extent = (y_max - y_min) * 1.2

    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_lines = lines1 + lines2

    rib_z_ranges = []
    for start, end in all_lines:
        z_range = get_rib_boundary_z_range(wing_mesh, start, end, params.plane_normal)
        if z_range is None:
            rib_z_ranges.append((None, None))
        else:
            rib_z_ranges.append(z_range)

    rib_union = create_holed_rib_union(wing_mesh, all_lines, params, params.plane_normal,
                                       params.rib_width, y_extent, params.cutout_margin,
                                       debug_dir=debug_dir, debug_format='step')
    if rib_union is None:
        print("No ribs generated – exiting.")
        return

    print("Subtracting ribs from wing...")
    # Only perform boolean if both are watertight
    if rib_union.is_volume and wing_mesh.is_watertight:
        try:
            cut_wing = trimesh.boolean.difference([wing_mesh, rib_union], engine='manifold')
        except Exception as e:
            print(f"Boolean difference failed: {e}")
            try:
                cut_wing = trimesh.boolean.difference([wing_mesh, rib_union], engine='scad')
            except:
                print("Scad fallback also failed; keeping original wing.")
                cut_wing = wing_mesh
    else:
        print("Skipping boolean difference: rib union or wing mesh is not watertight.")
        cut_wing = wing_mesh

    data_by_rib = collect_midpoints(all_lines, rib_z_ranges, params.z_step)
    bridge_meshes = create_bridge_meshes(data_by_rib, params.plane_normal,
                                         params.rib_width, params.bridge_height, y_extent)

    if bridge_meshes:
        print("Adding bridges...")
        bridge_union = None
        for b in bridge_meshes:
            if b.is_volume:
                if bridge_union is None:
                    bridge_union = b
                else:
                    try:
                        bridge_union = trimesh.boolean.union([bridge_union, b], engine='manifold')
                    except:
                        try:
                            bridge_union = trimesh.boolean.union([bridge_union, b], engine='scad')
                        except:
                            print("  Skipping a bridge due to union failure")
                            continue
        if bridge_union is not None:
            try:
                final_mesh = trimesh.boolean.union([cut_wing, bridge_union], engine='manifold')
            except Exception as e:
                print(f"Union failed: {e}")
                try:
                    final_mesh = trimesh.boolean.union([cut_wing, bridge_union], engine='scad')
                except:
                    print("Scad fallback failed; keeping wing without bridges.")
                    final_mesh = cut_wing
        else:
            final_mesh = cut_wing
    else:
        final_mesh = cut_wing

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
        input_stl_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin.stl",
        output_stl_path=r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin_out.stl",
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