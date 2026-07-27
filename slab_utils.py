# slab_utils.py – Mesh-optimised using trimesh (fully corrected)

import FreeCAD
import Part
import Mesh
import MeshPart
import math
import sys
import time
import os
import tempfile
import numpy as np
import trimesh
import shapely.geometry as sg
from shapely.ops import polygonize
from functools import wraps

# ------------------------------------------------------------
# Profiling (complete)
# ------------------------------------------------------------
class FreeCADProfiler:
    def __init__(self):
        self.operation_stack = []
        self.enabled = True

    def log_op(self, name: str, level: int = 0):
        if not self.enabled:
            return
        if self.operation_stack:
            prev_name, prev_start = self.operation_stack[-1]
            duration = time.time() - prev_start
            prefix = "  " * (level - 1)
            sys.stdout.write(f"{prefix}--> {prev_name} completed in {duration:.2f}s\n")
            sys.stdout.flush()
        self.operation_stack.append((name, time.time()))
        prefix = "  " * level
        sys.stdout.write(f"{prefix}Starting: {name}\n")
        sys.stdout.flush()

    def end_op(self, obj_count: int = 0):
        if not self.enabled or not self.operation_stack:
            return
        name, start_time = self.operation_stack.pop()
        duration = time.time() - start_time
        prefix = "  " * len(self.operation_stack)
        if obj_count > 0:
            sys.stdout.write(f"{prefix}✅ {name} completed in {duration:.2f}s ({obj_count} objects)\n")
        else:
            sys.stdout.write(f"{prefix}✅ {name} completed in {duration:.2f}s\n")
        sys.stdout.flush()

    def log(self, msg: str, level: int = 0):
        if not self.enabled:
            return
        prefix = "  " * level
        sys.stdout.write(f"{prefix}📝 {msg}\n")
        sys.stdout.flush()

profiler = FreeCADProfiler()
log = profiler.log
start_op = profiler.log_op
end_op = profiler.end_op
log_op = profiler.log

# ------------------------------------------------------------
# Mesh repair (from your pure Python script)
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# Convert FreeCAD shape <-> trimesh
# ------------------------------------------------------------
def shape_to_trimesh(shape, deflection=0.5, angular=0.5):
    """Convert a FreeCAD Part.Shape to a repaired trimesh.Trimesh."""
    if shape is None or shape.isNull():
        return None
    with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp:
        tmp_stl = tmp.name
    try:
        mesh_obj = MeshPart.meshFromShape(shape, LinearDeflection=deflection,
                                          AngularDeflection=angular, Relative=False)
        mesh_obj.write(tmp_stl)
        tm = trimesh.load(tmp_stl, force='mesh')
        os.unlink(tmp_stl)
        if tm is not None:
            tm = repair_mesh(tm)
        return tm
    except Exception as e:
        print(f"  shape_to_trimesh failed: {e}")
        try:
            os.unlink(tmp_stl)
        except:
            pass
        return None

def trimesh_to_freecad(mesh):
    """Convert a trimesh.Trimesh to FreeCAD Mesh.Mesh."""
    if mesh is None or len(mesh.vertices) == 0:
        return None
    with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp:
        tmp_stl = tmp.name
    mesh.export(tmp_stl)
    fc_mesh = Mesh.Mesh(tmp_stl)
    os.unlink(tmp_stl)
    return fc_mesh

# ------------------------------------------------------------
# Fast grid generation – all numpy (fixed TypeError)
# ------------------------------------------------------------
def create_angled_grid_lines(bb, params):
    """
    Returns two lists of (start, end) tuples, each a numpy array of shape (3,).
    """
    # Convert FreeCAD vectors to numpy arrays
    n = np.array(params.plane_normal)
    au = np.array(params.plane_axis_u)
    av = np.array(params.plane_axis_v)

    corners = [
        np.array([bb.XMin, bb.YMin, bb.ZMin]),
        np.array([bb.XMax, bb.YMin, bb.ZMin]),
        np.array([bb.XMin, bb.YMax, bb.ZMin]),
        np.array([bb.XMax, bb.YMax, bb.ZMin]),
        np.array([bb.XMin, bb.YMin, bb.ZMax]),
        np.array([bb.XMax, bb.YMin, bb.ZMax]),
        np.array([bb.XMin, bb.YMax, bb.ZMax]),
        np.array([bb.XMax, bb.YMax, bb.ZMax]),
    ]

    def dot(v, a): return np.dot(v, a)
    u_vals = [dot(c, au) for c in corners]
    v_vals = [dot(c, av) for c in corners]
    span_u = max(u_vals) - min(u_vals)
    span_v = max(v_vals) - min(v_vals)

    pd = params.primary_dir
    if pd is not None:
        pd = np.array(pd)
    else:
        pd = au if span_u >= span_v else av

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

    centre = np.array([(bb.XMin+bb.XMax)/2, (bb.YMin+bb.YMax)/2, (bb.ZMin+bb.ZMax)/2])
    line_len = np.linalg.norm([bb.XLength, bb.YLength, bb.ZLength]) * 2

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

# ------------------------------------------------------------
# Rib and bridge creation (from your pure Python script)
# ------------------------------------------------------------
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

def get_rib_boundary_z_range(wing_mesh, start, end, plane_normal, tol=1e-3):
    d = end - start
    length = np.linalg.norm(d)
    if length < 1e-6:
        return None
    d = d / length
    z_min = None
    z_max = None
    pts = [start + d * t * length for t in np.linspace(0, 1, 101)]
    contains = wing_mesh.contains(pts)
    for i, pt in enumerate(pts):
        if contains[i]:
            if z_min is None or pt[2] < z_min:
                z_min = pt[2]
            if z_max is None or pt[2] > z_max:
                z_max = pt[2]
    if z_min is not None and z_max is not None:
        return z_min, z_max
    return None

def create_holed_rib_union(wing_mesh, rib_lines, params, plane_normal, rib_width, y_extent, margin, debug_dir=None):
    log("  Building rib meshes...")
    rib_meshes = []
    for idx, (start, end) in enumerate(rib_lines):
        m = create_rib_mesh(start, end, plane_normal, rib_width, y_extent)
        if m is not None:
            rib_meshes.append(m)
    if not rib_meshes:
        return None
    log(f"  Unioning {len(rib_meshes)} ribs...")
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
    if not params.create_holes:
        return rib_union

    # ---- Cutouts ----
    log("  Generating cutouts...")
    wing_simple = wing_mesh.simplify_quadric_decimation(face_count=20000)
    cutout_prisms = []
    for idx, (start, end) in enumerate(rib_lines):
        d = end - start
        length = np.linalg.norm(d)
        if length < 1e-8:
            continue
        d = d / length
        n_rib = np.cross(d, plane_normal)
        n_rib = n_rib / np.linalg.norm(n_rib)
        origin = start

        try:
            path = trimesh.intersections.slice_mesh_plane(wing_simple, plane_normal=n_rib, plane_origin=origin)
        except:
            continue
        if path is None or len(path.vertices) < 3:
            continue

        verts = path.vertices
        edges = path.edges
        u = verts[edges[0][1]] - verts[edges[0][0]]
        u_norm = np.linalg.norm(u)
        if u_norm < 1e-8:
            u = np.cross(n_rib, [1,0,0])
            if np.linalg.norm(u) < 1e-8:
                u = np.cross(n_rib, [0,1,0])
            u = u / np.linalg.norm(u)
        else:
            u = u / u_norm
        v = np.cross(n_rib, u)
        v = v / np.linalg.norm(v)

        coords_2d = []
        for p in verts:
            dx = np.dot(p - origin, u)
            dy = np.dot(p - origin, v)
            coords_2d.append((dx, dy))

        lines_2d = []
        for edge in edges:
            p1 = coords_2d[edge[0]]
            p2 = coords_2d[edge[1]]
            lines_2d.append(sg.LineString([p1, p2]))
        merged = sg.MultiLineString(lines_2d)
        polygons = list(polygonize(merged))
        if not polygons:
            continue
        poly = max(polygons, key=lambda p: p.area)
        if poly.is_empty or poly.area < 1e-8:
            continue
        if not poly.is_valid:
            poly = poly.buffer(0)

        height = rib_width * 1.05
        transform = np.eye(4)
        transform[:3,:3] = np.column_stack([u, v, n_rib])
        transform[:3,3] = origin

        try:
            extruded = trimesh.creation.extrude_polygon(poly, height, transform=transform)
        except:
            continue
        if extruded is None or not extruded.is_volume:
            continue
        cutout_prisms.append(extruded)

    if cutout_prisms:
        log(f"  Unioning {len(cutout_prisms)} cutouts...")
        cutout_union = cutout_prisms[0]
        for cp in cutout_prisms[1:]:
            try:
                cutout_union = trimesh.boolean.union([cutout_union, cp], engine='manifold')
            except:
                try:
                    cutout_union = trimesh.boolean.union([cutout_union, cp], engine='scad')
                except:
                    print("  Skipping a cutout due to union failure")
                    continue
        if cutout_union.is_volume and rib_union.is_volume:
            log("  Subtracting cutouts from ribs...")
            try:
                rib_union = trimesh.boolean.difference([rib_union, cutout_union], engine='manifold')
            except:
                try:
                    rib_union = trimesh.boolean.difference([rib_union, cutout_union], engine='scad')
                except:
                    print("  Cutout subtraction failed – skipping holes.")
    return rib_union

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

def create_bridge_meshes(data_by_rib, plane_normal, rib_width, bridge_height):
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
            box = trimesh.creation.box(extents=[length, rib_width, bridge_height])
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

# ------------------------------------------------------------
# Main fast pipeline (manual profiling)
# ------------------------------------------------------------
def create_cut_result_mesh(wing_shape, params, doc, vis_cut_wing=False, vis_centre_lines=False, cutout_faces=None):
    start_op("create_cut_result_mesh")
    # Convert wing shape to trimesh
    wing_mesh = shape_to_trimesh(wing_shape, deflection=0.5, angular=0.5)
    if wing_mesh is None:
        raise RuntimeError("Failed to convert wing to trimesh")

    bb = wing_shape.BoundBox

    # 1. Grid lines (numpy arrays)
    lines1, lines2 = create_angled_grid_lines(bb, params)
    all_lines = lines1 + lines2  # list of (start, end) tuples

    # Convert to FreeCAD lines for visualisation and later use
    all_center_lines = []
    for start, end in all_lines:
        line = Part.makeLine(FreeCAD.Vector(*start), FreeCAD.Vector(*end))
        all_center_lines.append(line)

    # 2. Rib union with holes
    y_extent = (bb.YMax - bb.YMin) * 1.2
    rib_union = create_holed_rib_union(wing_mesh, all_lines, params,
                                       params.plane_normal, params.rib_width,
                                       y_extent, params.cutout_marging)
    if rib_union is None:
        print("No ribs generated – returning original wing.")
        end_op()
        return trimesh_to_freecad(wing_mesh), all_center_lines, []

    # 3. Subtract ribs from wing
    if wing_mesh.is_watertight and rib_union.is_volume:
        try:
            cut_wing = trimesh.boolean.difference([wing_mesh, rib_union], engine='manifold')
        except:
            try:
                cut_wing = trimesh.boolean.difference([wing_mesh, rib_union], engine='scad')
            except:
                cut_wing = wing_mesh
    else:
        cut_wing = wing_mesh

    # 4. Convert to FreeCAD mesh for visualisation
    cut_wing_fc = trimesh_to_freecad(cut_wing)
    if vis_cut_wing:
        from viz_utils import show_mesh
        show_mesh(cut_wing_fc, doc, "CutWingMesh", transparency=80)

    if vis_centre_lines:
        from viz_utils import show_rib_centre_lines
        show_rib_centre_lines(all_center_lines, doc)

    end_op()
    return cut_wing_fc, all_center_lines, [rib_union]

def create_bridges_trimmed_to_wing_mesh(data_by_rib, rib_center_lines, plane_normal,
                                        rib_width, bridge_height, wing_mesh,
                                        extend_length=5.0, doc=None, vis=False):
    start_op("create_bridges_trimmed_to_wing_mesh")
    # Convert FreeCAD lines to numpy points
    rib_lines = []
    for line in rib_center_lines:
        start = np.array(line.Vertexes[0].Point)
        end = np.array(line.Vertexes[-1].Point)
        rib_lines.append((start, end))

    # Convert data_by_rib (from point_utils) to our format: (z, point)
    bridge_data = {}
    for idx, data in data_by_rib.items():
        pts = [np.array(p) for p in data['mid']]
        if len(pts) < 2:
            continue
        pts_sorted = sorted(pts, key=lambda p: p[2])
        if len(pts_sorted) >= 2:
            dir0 = pts_sorted[1] - pts_sorted[0]
            dir0 = dir0 / np.linalg.norm(dir0)
            p_before = pts_sorted[0] - dir0 * extend_length
            dir_last = pts_sorted[-1] - pts_sorted[-2]
            dir_last = dir_last / np.linalg.norm(dir_last)
            p_after = pts_sorted[-1] + dir_last * extend_length
            pts_extended = [p_before] + pts_sorted + [p_after]
        else:
            continue
        bridge_data[idx] = [(p[2], p) for p in pts_extended]

    bridge_meshes = create_bridge_meshes(bridge_data, plane_normal, rib_width, bridge_height)
    if not bridge_meshes:
        end_op()
        return None

    # Union all bridges
    bridge_union = bridge_meshes[0]
    for bm in bridge_meshes[1:]:
        try:
            bridge_union = trimesh.boolean.union([bridge_union, bm], engine='manifold')
        except:
            try:
                bridge_union = trimesh.boolean.union([bridge_union, bm], engine='scad')
            except:
                print("  Skipping a bridge union due to failure")
                continue

    # Intersect with wing mesh (if wing_mesh is provided as trimesh)
    if wing_mesh is not None:
        try:
            trimmed = bridge_union.intersection(wing_mesh, engine='manifold')
        except:
            try:
                trimmed = bridge_union.intersection(wing_mesh, engine='scad')
            except:
                trimmed = bridge_union
        bridge_union = trimmed

    if vis:
        from viz_utils import show_mesh
        fc_mesh = trimesh_to_freecad(bridge_union)
        if fc_mesh:
            show_mesh(fc_mesh, doc, "BridgesMesh", color=(0.0,0.8,0.0), transparency=30)

    end_op()
    return trimesh_to_freecad(bridge_union)

def merge_and_show_final_mesh(cut_wing_fc, bridges_fc, doc, vis=True):
    start_op("merge_and_show_final_mesh")
    if bridges_fc is None or bridges_fc.CountFacets == 0:
        final_fc = cut_wing_fc
    else:
        with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp1:
            tmp1_stl = tmp1.name
        with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp2:
            tmp2_stl = tmp2.name
        try:
            cut_wing_fc.write(tmp1_stl)
            bridges_fc.write(tmp2_stl)
            tm1 = trimesh.load(tmp1_stl, force='mesh')
            tm2 = trimesh.load(tmp2_stl, force='mesh')
            if tm1 and tm2:
                try:
                    final_tm = trimesh.boolean.union([tm1, tm2], engine='manifold')
                except:
                    try:
                        final_tm = trimesh.boolean.union([tm1, tm2], engine='scad')
                    except:
                        final_tm = tm1
                final_fc = trimesh_to_freecad(final_tm)
            else:
                final_fc = cut_wing_fc
        except Exception as e:
            print(f"  Final union failed: {e}")
            final_fc = cut_wing_fc
        finally:
            try:
                os.unlink(tmp1_stl)
                os.unlink(tmp2_stl)
            except:
                pass
    if vis:
        from viz_utils import show_mesh
        show_mesh(final_fc, doc, "FinalMesh", color=(0.8,0.8,0.8), transparency=50)
    end_op()
    return final_fc