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
from shapely.geometry import Polygon as ShapelyPolygon, LineString
import trimesh.path.polygons
from shapely.ops import split
from shapely.geometry import Polygon as ShapelyPolygon, LineString
from viz_utils import show_mesh
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

def create_rib_surfaces_trimesh(wing_mesh, rib_center_lines_np, plane_normal_np, doc=None, vis=False):
    """
    Create flat rectangular surfaces for each rib line, spanning the wing's bounding box.
    No slicing – just pure visualization of the rib planes.
    """
    surfaces = []
    print(f"Creating {len(rib_center_lines_np)} rib planes...")

    # Get wing bounding box
    bbox = wing_mesh.bounds  # [[minx, miny, minz], [maxx, maxy, maxz]]
    min_pt = bbox[0]
    max_pt = bbox[1]
    # Approximate size for the rectangle (use diagonal of bbox)
    diag = np.linalg.norm(max_pt - min_pt)

    for i, (start, end) in enumerate(rib_center_lines_np):
        origin = (start + end) / 2.0
        rib_dir = end - start
        length = np.linalg.norm(rib_dir)
        if length < 1e-8:
            continue
        rib_dir /= length

        # Plane normal: perpendicular to rib_dir and construction normal
        plane_normal = np.cross(rib_dir, plane_normal_np)
        norm = np.linalg.norm(plane_normal)
        if norm < 1e-8:
            continue
        plane_normal /= norm

        # Build an orthonormal basis for the plane
        ref = np.array([1., 0., 0.])
        if abs(plane_normal[0]) > 0.9:
            ref = np.array([0., 1., 0.])
        u = np.cross(plane_normal, ref)
        u /= np.linalg.norm(u)
        v = np.cross(plane_normal, u)

        # Create a rectangular patch (size = diag * 0.8 to fit within bbox)
        half = diag * 0.5
        corners = [
            origin - half * u - half * v,
            origin + half * u - half * v,
            origin + half * u + half * v,
            origin - half * u + half * v
        ]
        # Build a simple quad mesh (2 triangles)
        verts = np.array(corners)
        faces = np.array([[0, 1, 2], [0, 2, 3]])
        mesh = trimesh.Trimesh(vertices=verts, faces=faces)
        surfaces.append(mesh)
        print(f"Rib {i}: rectangle created")

    print(f"Total surfaces created: {len(surfaces)}")

    # ---- Visualise if requested ----
    if vis and doc and surfaces:
        try:
            from viz_utils import show_mesh
            from slab_utils import trimesh_to_freecad
            combined = trimesh.util.concatenate(surfaces)
            fc_mesh = trimesh_to_freecad(combined)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "RibSurfaces", color=(0.2, 0.6, 1.0), transparency=40)
                doc.recompute()
                print("Rib surfaces visualised in FreeCAD.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return surfaces

def clip_surfaces_to_solid(rib_center_lines_np, wing_mesh, plane_normal_np,
                           doc=None, vis=False):
    """
    Create the exact rib cross-section surfaces by intersecting the wing mesh
    with each rib plane.

    Returns
    -------
    clipped : list[trimesh.Trimesh]
        Planar triangulated meshes representing the clipped rib surfaces.
    """

    clipped = []
    valid_indices = []
    print(f"Clipping {len(rib_center_lines_np)} rib planes against wing...")

    for i, (start, end) in enumerate(rib_center_lines_np):

        # ------------------------------------------------------------------
        # Rib plane
        # ------------------------------------------------------------------
        origin = (start + end) * 0.5

        rib_dir = end - start
        L = np.linalg.norm(rib_dir)
        if L < 1e-8:
            continue
        rib_dir /= L

        rib_plane_normal = np.cross(rib_dir, plane_normal_np)
        nL = np.linalg.norm(rib_plane_normal)
        if nL < 1e-8:
            continue
        rib_plane_normal /= nL

        # ------------------------------------------------------------------
        # Exact section curve (no cutting of the mesh)
        # ------------------------------------------------------------------
        try:
            section = wing_mesh.section(
                plane_origin=origin,
                plane_normal=rib_plane_normal
            )
        except Exception as e:
            print(f"Plane {i}: section failed: {e}")
            continue

        if section is None:
            print(f"Plane {i}: no intersection")
            continue

        # ------------------------------------------------------------------
        # Convert to planar coordinates
        # ------------------------------------------------------------------
        try:
            planar, to_3d = section.to_planar()
        except Exception as e:
            print(f"Plane {i}: to_planar failed: {e}")
            continue

        # planar is a Path2D
        polygons = planar.polygons_full

        if not polygons:
            print(f"Plane {i}: no closed polygons")
            continue

        rib_meshes = []

        # ------------------------------------------------------------------
        # Triangulate each polygon and transform back to 3D
        # ------------------------------------------------------------------
        for poly in polygons:

            if poly.is_empty or poly.area < 1e-8:
                continue

            try:
                verts2d, faces = trimesh.creation.triangulate_polygon(poly)
            except Exception as e:
                print(f"Plane {i}: triangulation failed: {e}")
                continue

            if len(verts2d) == 0 or len(faces) == 0:
                continue

            # Homogeneous 2D -> 3D transform
            verts_h = np.column_stack([
                verts2d,
                np.zeros(len(verts2d)),
                np.ones(len(verts2d))
            ])

            verts3d = (to_3d @ verts_h.T).T[:, :3]

            mesh = trimesh.Trimesh(
                vertices=verts3d,
                faces=faces,
                process=False
            )

            mesh.fix_normals()
            rib_meshes.append(mesh)

        if not rib_meshes:
            continue

        rib_mesh = trimesh.util.concatenate(rib_meshes)
        clipped.append(rib_mesh)
        valid_indices.append(i)

        print(f"Plane {i}: clipped surface created ({len(rib_mesh.vertices)} verts)")

    print(f"Total clipped surfaces created: {len(clipped)}")

    # ----------------------------------------------------------------------
    # Visualise
    # ----------------------------------------------------------------------
    if vis and doc and clipped:

        try:
            from viz_utils import show_mesh

            combined = trimesh.util.concatenate(clipped)
            fc_mesh = trimesh_to_freecad(combined)

            if fc_mesh:
                show_mesh(
                    fc_mesh,
                    doc,
                    "ClippedRibSurfaces",
                    color=(0.8, 0.2, 0.2),
                    transparency=20
                )

                doc.recompute()
                print("Clipped rib surfaces visualised in FreeCAD.")

        except Exception as e:
            print(f"Visualisation error: {e}")

    return clipped, valid_indices

def build_rib_segments_analytical(wing_mesh, all_lines_np, plane_normal_np,
                                  doc=None, vis=False):
    """
    Analitically build rib segments by intersecting a pre-built grid of
    rectangles with the wing's exact cross‑section polygon.

    Steps:
    1. For each rib centre line, compute its intersection points with every
       other line → segment boundaries.
    2. Obtain the wing's cross‑section polygon on the rib plane using
       wing_mesh.section() and to_planar().
    3. Map all crossing points into the polygon's local 2D coordinates
       (using the transformation matrix from to_planar).
    4. For each interval between consecutive u-values, build a rectangle
       in that same 2D space (interval horizontally, full polygon v-range).
    5. Intersect rectangle with wing polygon, triangulate, and transform
       back to 3D via the original 4x4 matrix.
    """
    # Precompute line data
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
        """3D line‑line intersection, returns 3D point or None"""
        A = np.column_stack((d1, -d2))
        b = p2 - p1
        try:
            st, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)
        except np.linalg.LinAlgError:
            return None
        if rank < 2:
            return None
        s, t = st[0], st[1]
        p_intersect = p1 + s * d1
        if np.linalg.norm(p2 + t * d2 - p_intersect) > 1e-6:
            return None
        return p_intersect

    all_segments = []

    for i, (line_i, ld_i) in enumerate(zip(all_lines_np, lines_data)):
        if ld_i is None:
            continue

        p_i = ld_i['point']      # start of line i
        d_i = ld_i['dir']        # unit direction

        # ---- 1. Collect crossing parameters (s) along the line ----
        crossing_params = []
        for j, ld_j in enumerate(lines_data):
            if j == i or ld_j is None:
                continue
            if abs(np.dot(d_i, ld_j['dir'])) > 0.9999:
                continue   # parallel → skip
            pt = intersect_lines(p_i, d_i, ld_j['point'], ld_j['dir'])
            if pt is None:
                continue
            s = np.dot(pt - p_i, d_i)
            crossing_params.append(s)

        # Add bounding box extremes to catch the whole wing
        bbox = wing_mesh.bounds
        bbox_center = (bbox[0] + bbox[1]) / 2
        proj = [np.dot(np.array(c), d_i) for c in [bbox[0], bbox[1], bbox_center]]
        s_min = min(proj) - 10.0
        s_max = max(proj) + 10.0
        crossing_params.append(s_min)
        crossing_params.append(s_max)

        crossing_params = sorted(set(crossing_params))

        # ---- 2. Wing cross‑section polygon + local 2D transform ----
        N_i = np.cross(d_i, plane_normal_np)
        nL = np.linalg.norm(N_i)
        if nL < 1e-8:
            continue
        N_i /= nL
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
        wing_poly = polygons[0]
        if wing_poly.is_empty or wing_poly.area < 1e-8:
            continue

        # Extract 2D‑to‑3D transformation details
        # to_3d is a 4x4 matrix; we need its inverse to map 3D → local 2D
        to_3d_mat = np.array(to_3d)   # shape (4,4)
        # The last column is the plane origin (homogeneous)
        origin_3d = to_3d_mat[:3, 3]
        # The first two columns are the orthonormal basis vectors (u, v)
        u_poly = to_3d_mat[:3, 0]
        v_poly = to_3d_mat[:3, 1]
        # Build inverse projection: point 3D → [u, v] = [ dot(P-O, u_poly), dot(P-O, v_poly) ]
        def project_to_local(pt3d):
            delta = pt3d - origin_3d
            return np.array([np.dot(delta, u_poly), np.dot(delta, v_poly)])

        # ---- 3. Map the wing polygon to local 2D (already done, it's wing_poly) ----
        # wing_poly is in local 2D coordinates; we can query its v‑range
        minx, miny, maxx, maxy = wing_poly.bounds
        v_min = miny
        v_max = maxy

        # ---- 4. Project crossing points into local 2D ----
        # crossing_params are s values along d_i. We must map 3D points
        # p_i + s * d_i to local 2D to get their u‑coordinate.
        crossing_u = []
        for s in crossing_params:
            pt3d = p_i + s * d_i
            uv = project_to_local(pt3d)
            crossing_u.append(uv[0])   # u‑coordinate in polygon's space

        crossing_u = sorted(set(crossing_u))

        # ---- 5. For each interval, build a rectangle in polygon's 2D space ----
        for k in range(len(crossing_u)-1):
            u0 = crossing_u[k]
            u1 = crossing_u[k+1]
            if u1 - u0 < 1e-8:
                continue

            # rectangle: [u0,u1] x [v_min, v_max]
            rect_pts = [
                (u0, v_min),
                (u1, v_min),
                (u1, v_max),
                (u0, v_max)
            ]
            rect = ShapelyPolygon(rect_pts)

            # Intersect with wing polygon
            try:
                intersection = rect.intersection(wing_poly)
            except Exception:
                continue
            if intersection.is_empty:
                continue

            if intersection.geom_type == 'Polygon':
                cells = [intersection]
            elif intersection.geom_type == 'MultiPolygon':
                cells = list(intersection.geoms)
            elif intersection.geom_type == 'GeometryCollection':
                cells = [g for g in intersection.geoms if g.geom_type == 'Polygon']
            else:
                continue

            for cell in cells:
                if cell.is_empty or cell.area < 1e-8:
                    continue
                try:
                    verts2d, faces = trimesh.creation.triangulate_polygon(cell)
                except Exception:
                    continue
                if len(verts2d) == 0 or len(faces) == 0:
                    continue

                # Map back to 3D using the same to_3d matrix
                verts_h = np.column_stack([verts2d, np.zeros(len(verts2d)), np.ones(len(verts2d))])
                verts3d = (to_3d_mat @ verts_h.T).T[:, :3]

                seg_mesh = trimesh.Trimesh(vertices=verts3d, faces=faces, process=False)
                seg_mesh.fix_normals()
                all_segments.append(seg_mesh)

    print(f"Created {len(all_segments)} rib segments (analytical grid).")

    # Visualisation
    if vis and doc and all_segments:
        try:
            combined = trimesh.util.concatenate(all_segments)
            fc_mesh = trimesh_to_freecad(combined)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "RibSegmentsAnalytical",
                          color=(0.2, 0.8, 0.4), transparency=30)
                doc.recompute()
                print("Rib segments visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return all_segments
