import numpy as np
import trimesh
from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union


def _face_plane_key(normal, point, normal_tol=6, dist_tol=4):
    n = tuple(np.round(normal, normal_tol))
    d = round(float(np.dot(point, normal)), dist_tol)
    return (n, d)


def merge_coplanar_faces(mesh, min_group_size=4):
    """
    Collapse groups of coplanar triangles sharing a plane (e.g. the
    hundreds of tiny wall triangles from unioning many small extruded
    prisms) into far fewer triangles, by re-triangulating each planar
    patch's merged boundary instead of keeping every original triangle.

    Lossless for genuinely flat regions -- not an approximation like
    quadric decimation. Groups below min_group_size are left alone.
    """
    face_normals = mesh.face_normals
    face_centers = mesh.triangles_center

    keys = {}
    for i in range(len(mesh.faces)):
        key = _face_plane_key(face_normals[i], face_centers[i])
        keys.setdefault(key, []).append(i)

    keep_mask = np.ones(len(mesh.faces), dtype=bool)
    new_verts, new_faces = [], []
    vert_offset = len(mesh.vertices)
    merged_count = 0

    for key, face_ids in keys.items():
        if len(face_ids) < min_group_size:
            continue
        face_ids = np.array(face_ids)
        normal = np.array(key[0])
        nlen = np.linalg.norm(normal)
        if nlen < 1e-8:
            continue
        normal /= nlen

        ref = np.array([1.0, 0.0, 0.0])
        if abs(normal[0]) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        axis_u = np.cross(normal, ref); axis_u /= np.linalg.norm(axis_u)
        axis_v = np.cross(normal, axis_u)
        origin = mesh.vertices[mesh.faces[face_ids[0]][0]]

        polys = []
        for fid in face_ids:
            tri = mesh.vertices[mesh.faces[fid]]
            pts2d = [(np.dot(p - origin, axis_u), np.dot(p - origin, axis_v)) for p in tri]
            poly = ShapelyPolygon(pts2d)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if not poly.is_empty and poly.area > 1e-12:
                polys.append(poly)
        if not polys:
            continue

        merged = unary_union(polys).buffer(0)
        if merged.is_empty:
            continue
        sub_polys = [merged] if merged.geom_type == 'Polygon' else list(merged.geoms)

        group_verts, group_faces, ok = [], [], True
        for poly in sub_polys:
            if poly.is_empty or poly.area < 1e-10:
                continue
            try:
                v2d, f = trimesh.creation.triangulate_polygon(poly)
            except Exception:
                ok = False
                break
            if len(v2d) == 0 or len(f) == 0:
                ok = False
                break
            pts3d = origin + v2d[:, 0:1] * axis_u + v2d[:, 1:2] * axis_v
            group_verts.append(pts3d)
            group_faces.append(f)  # indices relative to THIS poly's own v2d only

        if not ok or not group_faces:
            continue

        keep_mask[face_ids] = False
        for v, f in zip(group_verts, group_faces):
            new_verts.append(v)
            new_faces.append(f + vert_offset)
            vert_offset += len(v)
        merged_count += len(face_ids)

    if not new_faces:
        return mesh

    all_verts = np.vstack([mesh.vertices] + new_verts)
    all_faces = np.vstack([mesh.faces[keep_mask]] + new_faces)

    result = trimesh.Trimesh(vertices=all_verts, faces=all_faces, process=False)
    result.merge_vertices()
    try:
        result.update_faces(result.nondegenerate_faces())
    except Exception:
        pass
    result.fix_normals()

    print(f"merge_coplanar_faces: {len(mesh.faces)} -> {len(result.faces)} faces "
          f"({merged_count} originals merged)")
    return result