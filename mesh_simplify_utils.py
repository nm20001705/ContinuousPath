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


def drop_sliver_components(mesh, min_volume=1e-6, verbose=True,
                           max_fragment_volume=0.0):
    """
    Remove zero-volume debris the boolean leaves along the slit edges.

    Cutting 0.1mm-thick slabs out of a wing produces a lot of very thin
    geometry, and the boolean sheds a handful of degenerate scraps: on
    this wing, 7 extra "components" totalling 46 faces, every one of them
    with |volume| == 0.00 (largest 10.55 x 0.09 x 0.00 mm -- a flat
    ribbon with no thickness at all). They are harmless to slicers but
    they make the exported STL report watertight=False, which hides real
    problems.

    Components are dropped on VOLUME, not on size or count, so a
    genuinely detached but solid piece (which this design can legitimately
    produce, since full-depth slabs can isolate a cell) is always kept.
    """
    if mesh is None:
        return mesh
    try:
        comps = mesh.split(only_watertight=False)
    except Exception as e:
        if verbose:
            print(f"  sliver cleanup skipped ({e})")
        return mesh

    if len(comps) <= 1:
        return mesh

    keep = [c for c in comps if abs(c.volume) >= min_volume]

    # Anything still detached after the zero-volume debris is gone is a
    # solid crumb the bridges failed to hold on to. In this design every
    # rib slit is meant to be bridged, so a finished part should come out
    # as ONE body -- a detached piece means something else severed it,
    # typically a gap already present in the source model (on this
    # project's wingR1a, a malformed aileron hinge slot at Z~30 left two
    # crumbs of 7.31 and 0.03 mm^3). They print as loose specks inside the
    # part, so they are worth naming even when they are not removed.
    if len(keep) > 1:
        ranked = sorted(keep, key=lambda c: -abs(c.volume))
        total = sum(abs(c.volume) for c in ranked)
        crumbs = [c for c in ranked[1:]
                  if abs(c.volume) < max_fragment_volume]
        if verbose:
            print(f"  NOTE: {len(ranked)} disconnected bodies in the result "
                  f"(largest {abs(ranked[0].volume):.1f} mm^3 of {total:.1f}):")
            for c in ranked[1:6]:
                lo, hi = c.bounds
                print(f"    {abs(c.volume):10.2f} mm^3 at "
                      f"X {lo[0]:.1f}..{hi[0]:.1f}  Y {lo[1]:.1f}..{hi[1]:.1f}  "
                      f"Z {lo[2]:.1f}..{hi[2]:.1f}")
        if crumbs and max_fragment_volume > 0:
            keep = [c for c in keep if c not in crumbs]
            if verbose:
                print(f"  dropping {len(crumbs)} detached crumb(s) under "
                      f"{max_fragment_volume} mm^3")

    dropped = len(comps) - len(keep)
    if not keep or dropped == 0:
        return mesh

    # Take the surviving components exactly as split() produced them. Do
    # NOT merge_vertices() here: around 0.1mm slits there are legitimately
    # distinct vertices a hair apart, and welding them creates
    # non-manifold edges -- that alone turned a watertight result into a
    # leaky one, even though every dropped component was a closed
    # zero-volume shell that removing could not possibly open.
    out = trimesh.util.concatenate(keep) if len(keep) > 1 else keep[0]

    # Dropping whole components must not open the parts we keep. Deleting
    # degenerate FACES would (a watertight mesh can need zero-area faces
    # to stay closed), which is why this only ever removes complete
    # components -- and verifies afterwards, reverting if it went wrong.
    if mesh.is_watertight and not out.is_watertight:
        if verbose:
            print("  sliver cleanup would break watertightness -- keeping "
                  "the mesh intact instead")
        return mesh

    if verbose:
        lost = len(mesh.faces) - len(out.faces)
        print(f"  dropped {dropped} zero-volume sliver component(s) "
              f"({lost} faces); {len(keep)} solid part(s) kept, "
              f"watertight={out.is_watertight}")
    return out


def strip_degenerate_faces(mesh, verbose=False):
    """
    Drop zero-area faces, but never at the cost of watertightness.

    Removing degenerate triangles looks like an unconditionally safe
    cleanup and is not. A tessellation can contain a sliver that carries
    no area yet is still the only thing joining its three neighbours;
    deleting it leaves a three-edge hole.

    That is not hypothetical: on two of this project's inputs exactly ONE
    such face existed, and stripping it flipped the mesh from watertight
    to leaky. That in turn dropped the run off the fast trimesh path onto
    the BREP one, where it hung -- while the untouched tessellation had
    been a perfectly good volume all along.

    So the cleanup is applied only when it does not open the mesh. If the
    input was already watertight, that is the best state available and
    nothing here can improve on it.
    """
    if mesh is None:
        return mesh
    try:
        candidate = mesh.copy()
        candidate.update_faces(candidate.nondegenerate_faces())
    except Exception:
        return mesh

    if mesh.is_watertight and not candidate.is_watertight:
        if verbose:
            n = len(mesh.faces) - len(candidate.faces)
            print(f"  keeping {n} degenerate face(s): removing them would "
                  f"open an otherwise watertight mesh")
        return mesh
    return candidate
