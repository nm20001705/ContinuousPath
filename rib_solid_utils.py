# rib_solid_utils.py – turns the flat rib_segment meshes from
# build_rib_segments_analytical into a single solid, using the same
# per-line (prim, line_dir_3d, slab_normal) frame that bridges/holes use.

import trimesh
from viz_utils import show_mesh
from slab_utils import trimesh_to_freecad
from rib_slice_core import basis_vectors, prepare_segment
from solidify_utils import solidify_flat_segment, tree_union


def build_rib_solid_analytical(rib_segment_meshes, bridge_segments, primary_dir_np,
                                thickness, doc=None, vis=False, boolean_engine='manifold'):
    """
    rib_segment_meshes / bridge_segments : index-for-index, same lists
        used everywhere else (from build_rib_segments_analytical /
        segment_bounds).
    thickness : wall thickness of the rib material itself.

    Returns rib_solid (trimesh.Trimesh) or None.
    """
    if len(rib_segment_meshes) != len(bridge_segments):
        raise ValueError(
            f"rib_segment_meshes ({len(rib_segment_meshes)}) and "
            f"bridge_segments ({len(bridge_segments)}) must be the same "
            f"length and correspond index-for-index."
        )

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