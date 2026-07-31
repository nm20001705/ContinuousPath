# bridge_utils_analytical.py – analytical bridge generation using
# rib_segment meshes directly, built on the same shared slicing/void
# logic as hole_utils_analytical.py (see rib_slice_core.py).
#
# Width policy: constant bridge_height, centered in the chosen piece; if
# it doesn't fit, no bridge at that slice (no clamping). See
# rib_slice_core.bridge_width_interval.
#
# If thickness is given, ALSO builds the real extruded volume (per rib
# line, then unioned) in the same pass.

import numpy as np
import trimesh
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad
from rib_slice_core import (
    basis_vectors, iter_solid_pieces, t_to_3d, ribbons_to_mesh,
    bridge_width_interval, collect_line_intervals,
)
from solidify_utils import solidify_rib_line, tree_union


def create_bridges_analytical(rib_segment_meshes, bridge_segments, primary_dir_np,
                               z_vals, bridge_height, margin=0.0,
                               doc=None, vis=False,
                               thickness=None, boolean_engine='manifold'):
    """
    Returns (bridge_mesh, bridge_solid):
        bridge_mesh  : the flat visualization ribbon (or None).
        bridge_solid : the unioned extruded volume (or None if thickness
                        was not given, or nothing was generated).
    """
    if len(rib_segment_meshes) != len(bridge_segments):
        raise ValueError(
            f"rib_segment_meshes ({len(rib_segment_meshes)}) and "
            f"bridge_segments ({len(bridge_segments)}) must be the same "
            f"length and correspond index-for-index."
        )

    prim, u_ax, v_ax = basis_vectors(primary_dir_np)

    all_vertices = []
    all_faces = []
    vert_offset = 0

    solids = [] if thickness is not None else None

    def policy(piece):
        return bridge_width_interval(piece, bridge_height, margin)

    for seg_mesh, seg in zip(rib_segment_meshes, bridge_segments):
        ribbons = []

        for piece in iter_solid_pieces(seg_mesh, seg, prim, u_ax, v_ax, z_vals):
            result = policy(piece)
            if result is None:
                continue
            bridge_start, bridge_end = result

            pt_left = t_to_3d(bridge_start, piece['P_uv'], piece['t_rib'],
                               piece['bridge_dir_2d'], piece['plane_origin'], u_ax, v_ax)
            pt_right = t_to_3d(bridge_end, piece['P_uv'], piece['t_rib'],
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
                    intervals, frame['plane_offset'], frame['prim'],
                    frame['line_dir_3d'], frame['slab_normal'], thickness
                )
                if solid is not None:
                    solids.append(solid)

    bridge_mesh = None
    if all_vertices:
        verts_arr = np.array(all_vertices)
        faces_arr = np.array(all_faces)
        bridge_mesh = trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False)
        bridge_mesh.merge_vertices()
        bridge_mesh.fix_normals()
        print(f"Bridge mesh (analytical): {len(bridge_mesh.vertices)} verts, {len(bridge_mesh.faces)} faces")
    else:
        print("No bridge ribbon geometry generated.")

    bridge_solid = None
    if thickness is not None:
        bridge_solid = tree_union(solids, engine=boolean_engine)
        if bridge_solid is not None:
            print(f"Bridge solid: {len(bridge_solid.vertices)} verts, {len(bridge_solid.faces)} faces, "
                  f"watertight={bridge_solid.is_watertight}")
        else:
            print("No bridge solid geometry generated.")

    if vis and doc and bridge_mesh is not None:
        try:
            fc_mesh = trimesh_to_freecad(bridge_mesh)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "Bridges",
                          color=(0.9, 0.7, 0.1), transparency=20)
                doc.recompute()
                print("Bridges visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")

    return bridge_mesh, bridge_solid