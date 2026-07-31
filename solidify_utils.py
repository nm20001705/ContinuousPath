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
from shapely.geometry import box, Polygon as ShapelyPolygon

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

    Builds ONE continuous strip polygon from the (d0,d1,t_lo,t_hi)
    intervals -- walking the lower boundary forward and the upper
    boundary back -- instead of unioning a pile of independent
    axis-aligned rectangles. The old box-union approach staircased at
    every z_step wherever the width tapered, which both looked wrong
    and fed lots of tiny near-duplicate edges into the boolean engine.
    """
    if not intervals:
        return None

    # Sort by d0 to guarantee a well-ordered boundary walk regardless of
    # caller order.
    intervals = sorted(intervals, key=lambda iv: iv[0])

    lower = [(d0, t_lo) for d0, d1, t_lo, t_hi in intervals]
    lower.append((intervals[-1][1], intervals[-1][2]))
    upper = [(d0, t_hi) for d0, d1, t_lo, t_hi in intervals]
    upper.append((intervals[-1][1], intervals[-1][3]))

    ring = lower + list(reversed(upper))
    footprint = ShapelyPolygon(ring)
    if not footprint.is_valid:
        footprint = footprint.buffer(0)  # cheap self-intersection cleanup
    if footprint.is_empty or footprint.area < 1e-9:
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
        transform[:3, 3] = origin_const - (thickness / 2.0 + over_extrude) * slab_normal

        try:
            solid = trimesh.creation.extrude_triangulation(
                vertices=verts2d, faces=faces,
                height=thickness + 2 * over_extrude, transform=transform
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

    # Defensive: collapse any coincident-but-index-distinct vertices
    # (e.g. around a cutout rim from triangulate_polygon) before
    # extrusion -- extrude_triangulation finds its boundary loop by
    # edge-occurrence count, which duplicate vertices silently break.
    flat = trimesh.Trimesh(
        vertices=np.column_stack([verts2d, np.zeros(len(verts2d))]),
        faces=seg_mesh.faces,
        process=False,
    )
    flat.merge_vertices()
    try:
        flat.update_faces(flat.nondegenerate_faces())
    except Exception:
        pass
    verts2d = flat.vertices[:, :2]
    faces = flat.faces

    transform = np.eye(4)
    transform[:3, 0] = axis_u
    transform[:3, 1] = axis_v
    transform[:3, 2] = slab_normal
    transform[:3, 3] = plane_offset * slab_normal - (thickness / 2.0) * slab_normal

    try:
        solid = trimesh.creation.extrude_triangulation(
            vertices=verts2d, faces=faces, height=thickness, transform=transform
        )
    except Exception as e:
        print(f"solidify_flat_segment: extrude_triangulation failed: {e}")
        return None

    try:
        solid.merge_vertices()
    except Exception:
        pass
    try:
        flat.update_faces(flat.nondegenerate_faces())
    except Exception:
        pass
    try:
        solid.fix_normals()
    except Exception:
        pass
    try:
        trimesh.repair.fill_holes(solid, max_hole=thickness * 4)
    except Exception:
        pass
    try:
        solid.fix_normals()
    except Exception:
        pass

    if not solid.is_watertight or not solid.is_winding_consistent:
        print(f"solidify_flat_segment: still not clean after repair "
              f"(watertight={solid.is_watertight}, winding={solid.is_winding_consistent}, "
              f"verts={len(solid.vertices)}, faces={len(solid.faces)}) -- discarding")
        return None
    if solid.volume < 0:
        solid.invert()

    return solid
