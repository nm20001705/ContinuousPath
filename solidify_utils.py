# solidify_utils.py
#
# Turns a flat rib-line ribbon (bridge or hole strip, currently zero
# thickness) into an actual solid by extruding perpendicular to the rib's
# own plane (along slab_normal). Also provides a tree-reduction boolean
# union, since unioning hundreds of small prisms sequentially is both
# slow and numerically worse-conditioned than pairwise reduction.

import numpy as np
import trimesh
from shapely.geometry import box
from shapely.ops import unary_union

def solidify_rib_line(intervals, origin_const, axis_d, line_dir_3d, slab_normal,
                       thickness, over_extrude=0.0):
    """
    intervals : list of (d0, d1, t_start, t_end) in this rib line's own
        (d, t) parametrization -- d indexes along the rib's actual
        centerline (via axis_d = dir_rib/denom), NOT along the global
        stacking axis prim. Using prim directly here was the bug:
        d*prim moves vertically in global Z regardless of how the rib
        line itself is tilted, so the extruded strip ended up in the
        wrong plane entirely for any non-axis-aligned rib line.
    origin_const, axis_d : together give 3D(d,t) = origin_const +
        d*axis_d + t*line_dir_3d -- the same parametrization already
        used correctly by iter_solid_pieces/t_to_3d for the flat ribbon.
    """
    if not intervals:
        return None

    rects = [box(d0, t0, d1, t1) for d0, d1, t0, t1 in intervals if d1 > d0 and t1 > t0]
    if not rects:
        return None

    footprint = unary_union(rects)
    if footprint.is_empty:
        return None

    polys = [footprint] if footprint.geom_type == 'Polygon' else list(footprint.geoms)

    solids = []
    for poly in polys:
        if poly.is_empty or poly.area < 1e-9:
            continue
        try:
            verts2d, faces = trimesh.creation.triangulate_polygon(poly)
        except Exception:
            continue
        if len(verts2d) == 0 or len(faces) == 0:
            continue

        transform = np.eye(4)
        transform[:3, 0] = axis_d
        transform[:3, 1] = line_dir_3d
        transform[:3, 2] = slab_normal
        transform[:3, 3] = origin_const - (thickness / 2.0) * slab_normal

        try:
            solid = trimesh.creation.extrude_triangulation(
                vertices=verts2d, faces=faces, height=thickness + 2 * over_extrude, transform=transform
            )
        except Exception:
            continue
        solid.fix_normals()
        solids.append(solid)

    if not solids:
        return None
    result = trimesh.util.concatenate(solids) if len(solids) > 1 else solids[0]

    if not result.is_watertight or not result.is_winding_consistent:
        return result  # let tree_union's own guard catch/report this
    if result.volume < 0:
        result.invert()
    return result

def tree_union(meshes, engine='manifold'):
    """
    Boolean-union a list of meshes via pairwise tree reduction instead of
    sequential folding (reduce(a, b) -> reduce(result, c) -> ...).

    Sequential folding does O(n) boolean calls, each one against a mesh
    that keeps growing -- slow, and numerically the worst case (a huge
    accumulated mesh booleaned against one tiny sliver, repeated
    hundreds of times). Tree reduction does O(log n) calls, each between
    two similarly-sized operands.

    engine='manifold' requires `pip install manifold3d` and is currently
    the most robust + fastest correct trimesh boolean backend for this
    kind of many-small-coplanar-solids scenario. Avoid 'blender' (slow,
    subprocess-based) and be cautious with 'scad' (known non-manifold
    edge cases) for this use case.
    """
    meshes = [m for m in meshes if m is not None and len(m.vertices) > 0]
    if not meshes:
        return None
    while len(meshes) > 1:
        nxt = []
        for i in range(0, len(meshes) - 1, 2):
            nxt.append(trimesh.boolean.union([meshes[i], meshes[i + 1]], engine=engine))
        if len(meshes) % 2 == 1:
            nxt.append(meshes[-1])
        meshes = nxt
    return meshes[0]

def solidify_flat_segment(seg_mesh, dir_rib, slab_normal, plane_offset, thickness):
    """
    Extrude an already-planar rib segment mesh (from
    build_rib_segments_analytical) into a solid, thickness/2 on each
    side of its own plane along slab_normal.

    The segment lies in the plane spanned by dir_rib and the in-plane
    direction perpendicular to it (normal = slab_normal) -- NOT the
    (prim, line_dir_3d) frame used elsewhere for per-Z-slice bridge/hole
    bookkeeping. Reusing that frame here was the bug: prim (global Z)
    is generally not within this segment's own tilted plane, so
    reprojecting onto (prim, line_dir_3d) sheared every rib segment
    toward the same global direction, making them all look parallel
    instead of following their own rib line's angle.
    """
    if seg_mesh is None or len(seg_mesh.vertices) == 0:
        return None

    axis_u = dir_rib / np.linalg.norm(dir_rib)
    axis_v = np.cross(slab_normal, axis_u)
    nv = np.linalg.norm(axis_v)
    if nv < 1e-8:
        return None
    axis_v /= nv

    verts2d = np.column_stack([
        seg_mesh.vertices @ axis_u,
        seg_mesh.vertices @ axis_v,
    ])

    transform = np.eye(4)
    transform[:3, 0] = axis_u
    transform[:3, 1] = axis_v
    transform[:3, 2] = slab_normal
    transform[:3, 3] = plane_offset * slab_normal - (thickness / 2.0) * slab_normal

    try:
        solid = trimesh.creation.extrude_triangulation(
            vertices=verts2d, faces=seg_mesh.faces, height=thickness, transform=transform
        )
    except Exception:
        return None
    solid.fix_normals()

    if not solid.is_watertight or not solid.is_winding_consistent:
        return None
    if solid.volume < 0:
        solid.invert()

    return solid
