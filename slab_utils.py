# slab_utils.py – Mesh-optimised using trimesh (fully corrected)

import FreeCAD
import Part
import Mesh
import MeshPart
import math
import os
import tempfile
import numpy as np
import trimesh
from shapely.geometry import Polygon as ShapelyPolygon
import trimesh.path.polygons
from shapely.geometry import Polygon as ShapelyPolygon
from viz_utils import show_mesh
from shapely.ops import linemerge, polygonize
from shapely.geometry import Polygon as ShapelyPolygon, LineString
from shapely.affinity import affine_transform
from scipy.spatial import ConvexHull

def precompute_slices(wing_mesh, prim, z_step, d_min, d_max):
    """Return dict { rounded_d : (polygon_2d, to_3d_matrix) } for every slice."""
    # Same basis vectors that bridge/hole functions will use
    if abs(prim[0]) > 0.9:
        u_ax = np.cross(prim, [0, 1, 0])
    else:
        u_ax = np.cross(prim, [1, 0, 0])
    u_ax = u_ax / np.linalg.norm(u_ax)
    v_ax = np.cross(prim, u_ax)

    slices = {}
    d = d_min
    while d <= d_max + 1e-9:
        plane_origin = d * prim
        try:
            result = trimesh.intersections.mesh_plane(
                wing_mesh, plane_normal=prim, plane_origin=plane_origin
            )
        except Exception:
            d += z_step
            continue
        if isinstance(result, tuple):
            lines = result[0]
        else:
            lines = result
        if lines is None or len(lines) == 0:
            d += z_step
            continue

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
            d += z_step
            continue

        merged = linemerge(segments_2d)
        if merged.is_empty:
            d += z_step
            continue

        polys = list(polygonize(merged))
        if not polys:
            all_pts = []
            for line in segments_2d:
                all_pts.extend(line.coords)
            pts = np.array(all_pts)
            if len(pts) < 3:
                d += z_step
                continue
            try:
                hull = ConvexHull(pts)
                poly = ShapelyPolygon(pts[hull.vertices])
            except Exception:
                d += z_step
                continue
        else:
            poly = max(polys, key=lambda p: p.area)

        if poly.is_empty or poly.area < 1e-8:
            d += z_step
            continue

        to_3d = np.eye(4)
        to_3d[:3, 0] = u_ax
        to_3d[:3, 1] = v_ax
        to_3d[:3, 3] = origin

        slices[round(d, 6)] = (poly, to_3d)
        d += z_step

    return slices

# ------------------------------------------------------------
# Mesh repair (from your pure Python script)
# ------------------------------------------------------------
def repair_mesh(mesh):
    """
    Minimal repair that fixes only small defects while preserving
    intentional holes (servo cutouts, lightening holes).
    """
    if mesh is None:
        return None
    if isinstance(mesh, trimesh.Scene):
        meshes = [g for g in mesh.geometry.values() if isinstance(g, trimesh.Trimesh)]
        if not meshes:
            return None
        mesh = trimesh.util.concatenate(meshes) if len(meshes) > 1 else meshes[0]
        if mesh is None:
            return None
    if not isinstance(mesh, trimesh.Trimesh) or len(mesh.faces) == 0:
        return None

    # Safe, non‑shape‑altering operations
    try:
        mesh.merge_vertices()
    except:
        pass
    try:
        mesh.update_faces(mesh.nondegenerate_faces())
    except:
        pass
    try:
        mesh.fix_normals()
    except:
        pass

    # Fill only very small holes (≤0.1 mm) to fix mesh defects without closing
    # intentional cutouts which are typically much larger.
    try:
        mesh = trimesh.repair.fill_holes(mesh, max_hole=0.1)
    except:
        pass

    # Re‑apply normals after hole filling
    try:
        mesh.fix_normals()
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
    Returns
    -------
    all_segments : list[trimesh.Trimesh]
    segment_bounds : list[tuple]
        Each tuple is (line_index, s_start, s_end) where s is the parameter
        along the rib line defined by p0 + s * dir_rib.
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
    segment_bounds = []          # <-- new

    for i, (line_i, ld_i) in enumerate(zip(all_lines_np, lines_data)):
        if ld_i is None:
            continue

        p_i = ld_i['point']
        d_i = ld_i['dir']

        # ---- 1. Collect crossing parameters (s) along the line ----
        crossing_params = []
        for j, ld_j in enumerate(lines_data):
            if j == i or ld_j is None:
                continue
            if abs(np.dot(d_i, ld_j['dir'])) > 0.9999:
                continue
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
        wing_poly = max(polygons, key=lambda p: p.area)
        if wing_poly.is_empty or wing_poly.area < 1e-6:
            continue

        # section.to_planar() picks an ARBITRARY in-plane basis, but the
        # cell-splitting below cuts axis-aligned bands in u and therefore
        # only works if u runs along the rib. Whether it does is pure
        # luck: with construction plane XZ / primary_dir Z trimesh happens
        # to give u == d_i, but with XY / primary_dir Y it gives v == d_i
        # instead -- every crossing then projects to the SAME u, the band
        # widths collapse below the 1e-8 guard, and the rib yields zero
        # segments. So re-express the polygon in our own frame with u
        # pinned to the rib direction. As a bonus this makes u strictly
        # increasing in s, so crossing_u and crossing_params below stay
        # in lockstep instead of being two independently-sorted lists.
        to_3d_in = np.array(to_3d)      # (4,4) from to_planar
        origin_3d = to_3d_in[:3, 3]
        u_in = to_3d_in[:3, 0]
        v_in = to_3d_in[:3, 1]

        u_poly = d_i                    # in-plane, unit, along the rib
        v_poly = np.cross(N_i, u_poly)  # in-plane, unit (N_i ⟂ d_i)

        # to_planar's frame -> ours. Both share origin_3d and are
        # orthonormal, so this is a pure 2x2 rotation/reflection.
        a, b = np.dot(u_in, u_poly), np.dot(v_in, u_poly)
        c, e = np.dot(u_in, v_poly), np.dot(v_in, v_poly)
        wing_poly = affine_transform(wing_poly, [a, b, c, e, 0.0, 0.0])
        if wing_poly.is_empty or wing_poly.area < 1e-6:
            continue

        to_3d_mat = np.eye(4)
        to_3d_mat[:3, 0] = u_poly
        to_3d_mat[:3, 1] = v_poly
        to_3d_mat[:3, 2] = N_i
        to_3d_mat[:3, 3] = origin_3d

        def project_to_local(pt3d):
            delta = pt3d - origin_3d
            return np.array([np.dot(delta, u_poly), np.dot(delta, v_poly)])

        minx, miny, maxx, maxy = wing_poly.bounds
        v_min, v_max = miny, maxy

        # Map crossing parameters to local u coordinates
        crossing_u = []
        for s in crossing_params:
            pt3d = p_i + s * d_i
            uv = project_to_local(pt3d)
            crossing_u.append(uv[0])
        crossing_u = sorted(set(crossing_u))

        # ---- 3. For each interval, build a rectangle and intersect ----
        for k in range(len(crossing_u)-1):
            u0 = crossing_u[k]
            u1 = crossing_u[k+1]
            if u1 - u0 < 1e-8:
                continue

            rect = ShapelyPolygon([(u0, v_min), (u1, v_min),
                                   (u1, v_max), (u0, v_max)])
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

                verts_h = np.column_stack([verts2d, np.zeros(len(verts2d)), np.ones(len(verts2d))])
                verts3d = (to_3d_mat @ verts_h.T).T[:, :3]

                seg_mesh = trimesh.Trimesh(vertices=verts3d, faces=faces, process=False)
                seg_mesh.merge_vertices()   # collapses the duplicate rim vertices
                seg_mesh.fix_normals()
                all_segments.append(seg_mesh)   # <-- this line was missing

                # ----- Record the segment bounds -----
                segment_bounds.append((i, crossing_params[k], crossing_params[k+1]))

    print(f"Created {len(all_segments)} rib segments (analytical grid).")

    # Visualisation unchanged...
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

    return all_segments, segment_bounds


def merge_rib_segments_by_line(rib_segments, segment_bounds, lines_data):
    """
    Recombine the per-segment meshes from build_rib_segments_analytical
    back into one mesh per original rib line (segments are only split at
    rib-rib intersections; see the crossing_u logic above). Downstream
    slicing (iter_solid_pieces, used by create_bridges_analytical) then
    runs continuously across the full line instead of restarting at every
    intersection -- restarting is what was causing bridges to break right
    at rib crossings, since each segment's own Z-range stops just short of
    the shared boundary.

    Returns
    -------
    line_meshes : list[trimesh.Trimesh]
    line_segments : list[dict]
        One {'p0','p1','dir','slab_normal'} bridge-segment dict per line,
        index-aligned with line_meshes.
    """
    by_line = {}
    for seg_mesh, (line_idx, s0, s1) in zip(rib_segments, segment_bounds):
        by_line.setdefault(line_idx, []).append(seg_mesh)

    line_meshes = []
    line_segments = []
    for line_idx, meshes in by_line.items():
        ld = lines_data[line_idx]
        if ld is None:
            continue
        merged = trimesh.util.concatenate(meshes) if len(meshes) > 1 else meshes[0]
        merged.merge_vertices()
        merged.fix_normals()
        line_meshes.append(merged)
        line_segments.append({
            'p0': ld['point'],
            'p1': ld['point'] + ld['dir'],
            'dir': ld['dir'],
            'slab_normal': ld['slab_normal'],
        })

    return line_meshes, line_segments
