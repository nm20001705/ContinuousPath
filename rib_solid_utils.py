# rib_solid_utils.py -- builds the rib solid.
#
# build_full_rib_solid (new, recommended): builds each rib as ONE big
# rectangular slab, sized to comfortably clear the wing everywhere in
# its own plane, plus a fixed protrusion_margin -- instead of clipping
# to the wing's tessellated cross-section like the old segment-based
# approach did. Because every slab is guaranteed to extend past the
# real wing surface everywhere, the eventual wing_shape.cut(...) always
# cuts fully through the skin (no more coincident-face / exact-
# surface-match problems from the old approach, whose boundary
# conformed EXACTLY to the wing's own tessellation). It also sidesteps
# the whole class of cutout-interior-ring / duplicate-vertex bugs the
# segment-based version had, since these slabs are always simple,
# hole-free rectangles -- nothing to trip up extrude_triangulation's
# boundary-edge detection.
#
# build_rib_solid_from_segments (old): kept for reference/comparison --
# this is the previous per-segment approach.

import numpy as np
import trimesh
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad
from rib_slice_core import basis_vectors, prepare_segment
from solidify_utils import solidify_flat_segment, tree_union


def _bbox_corners(bbox):
    (xmin, ymin, zmin), (xmax, ymax, zmax) = bbox
    return np.array([
        [xmin, ymin, zmin], [xmin, ymin, zmax],
        [xmin, ymax, zmin], [xmin, ymax, zmax],
        [xmax, ymin, zmin], [xmax, ymin, zmax],
        [xmax, ymax, zmin], [xmax, ymax, zmax],
    ])


def build_full_rib_solid(all_lines_np, plane_normal_np, wing_mesh,
                          thickness, protrusion_margin=2.0,
                          doc=None, vis=False, boolean_engine='manifold'):
    """
    all_lines_np : the FULL rib center lines (e.g. lines1+lines2 from
        create_angled_grid_lines) -- NOT split at crossings, NOT
        clipped to the wing.
    protrusion_margin : extra distance (model units) added on every
        side beyond the wing's own bounding box, in BOTH in-plane
        directions of each rib's plane. Guarantees every slab edge
        clears the real wing surface by at least this much, everywhere
        -- the wing is entirely contained within its own bbox, so a
        rectangle sized to bbox+margin in both directions necessarily
        pokes out past the actual (tighter, curved) wing surface too.

    Returns rib_solid (trimesh.Trimesh) or None.
    """
    corners = _bbox_corners(wing_mesh.bounds)

    solids = []
    n_skipped = 0
    for start, end in all_lines_np:
        dir_rib = end - start
        length = np.linalg.norm(dir_rib)
        if length < 1e-8:
            n_skipped += 1
            continue
        dir_rib = dir_rib / length

        slab_normal = np.cross(dir_rib, plane_normal_np)
        nsn = np.linalg.norm(slab_normal)
        if nsn < 1e-8:
            n_skipped += 1
            continue
        slab_normal /= nsn

        # True orthonormal in-plane basis: dir_rib and axis_v are both
        # unit length, perpendicular to slab_normal AND to each other,
        # by construction. They don't need to align with any global
        # axis (prim, etc.) -- they just need to span the rib's plane,
        # which is all solidify_flat_segment needs downstream.
        axis_v = np.cross(slab_normal, dir_rib)

        mid = (start + end) / 2.0
        proj_u = np.dot(corners - mid, dir_rib)
        proj_v = np.dot(corners - mid, axis_v)

        u_min, u_max = proj_u.min() - protrusion_margin, proj_u.max() + protrusion_margin
        v_min, v_max = proj_v.min() - protrusion_margin, proj_v.max() + protrusion_margin
        u_center, half_u = (u_min + u_max) / 2.0, (u_max - u_min) / 2.0
        v_center, half_v = (v_min + v_max) / 2.0, (v_max - v_min) / 2.0

        # Moving along dir_rib and axis_v only (both perpendicular to
        # slab_normal) keeps `origin` exactly on the rib's original
        # plane -- so plane_offset below still lands at THIS specific
        # rib's correct physical location (preserving rib spacing).
        origin = mid + u_center * dir_rib + v_center * axis_v

        rect_corners = [
            origin - half_u * dir_rib - half_v * axis_v,
            origin + half_u * dir_rib - half_v * axis_v,
            origin + half_u * dir_rib + half_v * axis_v,
            origin - half_u * dir_rib + half_v * axis_v,
        ]
        verts = np.array(rect_corners)
        faces = np.array([[0, 1, 2], [0, 2, 3]])
        rect_mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

        plane_offset = np.dot(origin, slab_normal)
        solid = solidify_flat_segment(rect_mesh, dir_rib, slab_normal, plane_offset, thickness)
        if solid is None:
            print(f"[full_rib_solid] slab for line starting at {start} failed to solidify")
            n_skipped += 1
            continue
        solids.append(solid)

    print(f"[full_rib_solid] {len(solids)} rib slabs built, {n_skipped} lines skipped")
    rib_solid = tree_union(solids, engine=boolean_engine)
    if rib_solid is not None:
        print(f"Full rib solid: {len(rib_solid.vertices)} verts, {len(rib_solid.faces)} faces, "
              f"watertight={rib_solid.is_watertight}")
    else:
        print("No full rib solid geometry generated.")

    if vis and doc and rib_solid is not None:
        try:
            fc_mesh = trimesh_to_freecad(rib_solid)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "RibSolidFull", color=(0.4, 0.4, 0.9), transparency=20)
                doc.recompute()
                print("Full rib solid visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return rib_solid


def build_rib_solid_from_segments(rib_segment_meshes, bridge_segments, primary_dir_np,
                                   thickness, doc=None, vis=False, boolean_engine='manifold'):
    """
    OLD approach, kept for reference/comparison: builds the rib solid
    from the wing-clipped, crossing-split segments produced by
    build_rib_segments_analytical. Its outer boundary conforms EXACTLY
    to the wing's own tessellation, which can leave the eventual
    wing_shape.cut(...) unable to cut fully through the skin
    (coincident faces), and can silently drop segments near
    cutout-adjacent cells. Prefer build_full_rib_solid.
    """
    prim, u_ax, v_ax = basis_vectors(primary_dir_np)

    solids = []
    n_none = 0
    for idx, (seg_mesh, seg) in enumerate(zip(rib_segment_meshes, bridge_segments)):
        prep = prepare_segment(seg, prim, u_ax, v_ax)
        if prep is None:
            print(f"[rib_solid] segment {idx}: prepare_segment failed (degenerate direction)")
            n_none += 1
            continue
        solid = solidify_flat_segment(
            seg_mesh, prep['dir_rib'], prep['slab_normal'], prep['plane_offset'], thickness
        )
        if solid is None:
            print(f"[rib_solid] segment {idx}: solidify_flat_segment returned None "
                  f"(verts={len(seg_mesh.vertices)}, faces={len(seg_mesh.faces)})")
            n_none += 1
            continue
        solids.append(solid)

    print(f"[rib_solid] {len(solids)} solids built, {n_none} segments dropped")
    rib_solid = tree_union(solids, engine=boolean_engine)
    if rib_solid is not None:
        print(f"Rib solid: {len(rib_solid.vertices)} verts, {len(rib_solid.faces)} faces, "
              f"watertight={rib_solid.is_watertight}")
    else:
        print("No rib solid geometry generated.")

    if vis and doc and rib_solid is not None:
        try:
            fc_mesh = trimesh_to_freecad(rib_solid)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "RibSolid", color=(0.4, 0.4, 0.9), transparency=20)
                doc.recompute()
                print("Rib solid visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return rib_solid


# Backwards-compatible alias, in case anything still imports the old
# name -- now points explicitly at the segment-based implementation.
build_rib_solid_analytical = build_rib_solid_from_segments