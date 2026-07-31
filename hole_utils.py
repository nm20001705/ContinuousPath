# hole_utils_analytical.py – analytical hole generation using rib_segment
# meshes directly, built on the shared slicing/void logic in
# rib_slice_core.py.
#
# Width policy for holes: width = point_condition(x) * available_width,
# centered in the chosen solid piece, clamped never to exceed the
# available (margin-reduced) interval. See rib_slice_core.hole_width_interval.
#
# If thickness is given, ALSO builds the real extruded volume (per rib
# line, then unioned) in the same pass -- no need to re-walk every slice
# a second time just to solidify.

import numpy as np
import trimesh
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad
from rib_slice_core import (
    basis_vectors, iter_solid_pieces, t_to_3d, ribbons_to_mesh,
    hole_width_interval, collect_line_intervals,
)
from solidify_utils import solidify_rib_line, tree_union


def create_holes_analytical(rib_segment_meshes, bridge_segments, primary_dir_np,
                             z_vals, point_condition=None, hole_margin=0.0,
                             doc=None, vis=False,
                             thickness=None, boolean_engine='manifold',
                             over_extrude=0.02):
    """
    rib_segment_meshes / bridge_segments : see module docstring / earlier
        turns -- index i corresponds 1:1 (both from segment_bounds).
    z_vals : np.ndarray
        Global precomputed Z-slice positions.
    thickness : float or None
        If given, also extrude each rib line's hole strip into a real
        volume (perpendicular to that rib's own plane) and union them
        all into a single solid, returned as `hole_solid`. Also used to
        widen the effective hole_margin so the hole cutter stays clear
        of the crossing rib's real extruded volume (see
        rib_slice_core.hole_width_interval).

    Returns (hole_mesh, hole_solid):
        hole_mesh  : the flat visualization ribbon (or None).
        hole_solid : the unioned extruded volume (or None if thickness
                     was not given, or nothing was generated).
    """
    if point_condition is None:
        point_condition = lambda x: 1.0

    if len(rib_segment_meshes) != len(bridge_segments):
        raise ValueError(
            f"rib_segment_meshes ({len(rib_segment_meshes)}) and "
            f"bridge_segments ({len(bridge_segments)}) must be the same "
            f"length and correspond index-for-index."
        )

    # hole_width_interval needs a numeric thickness for its margin term
    # even when the caller didn't ask for a solid (thickness=None) --
    # fall back to 0 so the margin behaves exactly as before in that case.
    margin_thickness = thickness if thickness is not None else 0.0

    prim, u_ax, v_ax = basis_vectors(primary_dir_np)

    all_vertices = []
    all_faces = []
    vert_offset = 0

    solids = [] if thickness is not None else None

    def policy(piece):
        return hole_width_interval(piece, point_condition, hole_margin, margin_thickness)

    for seg_mesh, seg in zip(rib_segment_meshes, bridge_segments):
        ribbons = []

        for piece in iter_solid_pieces(seg_mesh, seg, prim, u_ax, v_ax, z_vals):
            result = policy(piece)
            if result is None:
                continue
            hole_start, hole_end = result

            pt_left = t_to_3d(hole_start, piece['P_uv'], piece['t_rib'],
                               piece['bridge_dir_2d'], piece['plane_origin'], u_ax, v_ax)
            pt_right = t_to_3d(hole_end, piece['P_uv'], piece['t_rib'],
                                piece['bridge_dir_2d'], piece['plane_origin'], u_ax, v_ax)

            ribbons.append((piece['d'], pt_left, pt_right))

        if len(ribbons) >= 2:
            ribbons.sort(key=lambda r: r[0])
            verts, faces = ribbons_to_mesh(ribbons)
            all_vertices.extend(verts)
            all_faces.extend([[a + vert_offset, b + vert_offset, c + vert_offset] for a, b, c in faces])
            vert_offset += len(verts)

        if thickness is not None:
            intervals, frame = collect_line_intervals(seg_mesh, seg, prim, u_ax, v_ax, z_vals, policy)
            if intervals and frame is not None:
                solid = solidify_rib_line(
                    intervals, frame['origin_const'], frame['axis_d'],
                    frame['line_dir_3d'], frame['slab_normal'], thickness,
                    over_extrude=over_extrude
                )
                if solid is not None:
                    solids.append(solid)

    hole_mesh = None
    if all_vertices:
        verts_arr = np.array(all_vertices)
        faces_arr = np.array(all_faces)
        hole_mesh = trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False)
        hole_mesh.merge_vertices()
        hole_mesh.fix_normals()
        print(f"Hole mesh (analytical): {len(hole_mesh.vertices)} verts, {len(hole_mesh.faces)} faces")
    else:
        print("No hole ribbon geometry generated.")

    hole_solid = None
    if thickness is not None:
        hole_solid = tree_union(solids, engine=boolean_engine)
        if hole_solid is not None:
            print(f"Hole solid: {len(hole_solid.vertices)} verts, {len(hole_solid.faces)} faces, "
                  f"watertight={hole_solid.is_watertight}")
        else:
            print("No hole solid geometry generated.")

    if vis and doc and hole_mesh is not None:
        try:
            fc_mesh = trimesh_to_freecad(hole_mesh)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "Holes",
                          color=(0.1, 0.1, 0.1), transparency=80)
                doc.recompute()
                print("Holes visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return hole_mesh, hole_solid