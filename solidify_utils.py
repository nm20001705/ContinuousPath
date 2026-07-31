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


def solidify_rib_line(intervals, plane_offset, prim, line_dir_3d, slab_normal, thickness):
    """
    intervals : list of (d0, d1, t_start, t_end)
        The rectangles in this rib line's own local (d, t) plane -- one
        per pair of adjacent Z-slices. This is the same data that used
        to go straight into a quad strip; keeping it in 2D first (rather
        than projecting to 3D immediately) is what lets us hand it to
        shapely for clean polygon unioning before triangulating.
    plane_offset : float
        dot(p0, slab_normal) for this rib line -- the constant position
        of this rib's plane along its own normal.
    prim, line_dir_3d, slab_normal : np.ndarray (3,)
        Orthonormal local frame for this rib line: prim is the Z-stacking
        direction, line_dir_3d is the in-plane "width" direction (this is
        the same vector bridge_dir_2d was derived from, kept here as a
        full 3D unit vector instead of projected into u_ax/v_ax).
    thickness : float
        Wall thickness, extruded symmetrically about the rib plane
        (+-thickness/2 along slab_normal).

    Returns a single trimesh.Trimesh (possibly multiple disjoint islands
    concatenated together) or None if there's no geometry.
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

        # Local frame -> world. x=d axis (prim), y=t axis (line_dir_3d),
        # z=extrusion axis (slab_normal). Centered on the rib plane.
        transform = np.eye(4)
        transform[:3, 0] = prim
        transform[:3, 1] = line_dir_3d
        transform[:3, 2] = slab_normal
        transform[:3, 3] = plane_offset * slab_normal - (thickness / 2.0) * slab_normal

        solid = trimesh.creation.extrude_triangulation(
            vertices=verts2d, faces=faces, height=thickness, transform=transform
        )
        solid.fix_normals()
        solids.append(solid)

    if not solids:
        return None
    return trimesh.util.concatenate(solids) if len(solids) > 1 else solids[0]


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