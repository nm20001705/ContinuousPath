# assembly_utils_freecad.py -- final wing assembly.
#
# Strategy:
#   1. rib_solid - bridge_solid - hole_solid  -> done ENTIRELY in trimesh.
#      All three are confirmed valid volumes (watertight, positive
#      volume), so this boolean is fast and doesn't need BREP at all.
#   2. Simplify the resulting single mesh (it's mostly flat extruded
#      prism walls -- very over-triangulated for its actual shape).
#   3. Convert ONLY that one simplified mesh to a Part solid via
#      makeShapeFromMesh + sew. This is the slow step, but now it only
#      runs once, on one (much smaller) mesh, instead of three times.
#   4. wing_shape.cut(cut_part) -- the one real BREP boolean, done
#      against the actual wing geometry because wing_shape has
#      topological defects that make it permanently invalid as a
#      trimesh volume (no non-destructive amount of hole-filling fixes
#      that). FreeCAD's own boolean tolerates it fine.

import trimesh
import Part
from slab_utils import trimesh_to_freecad
from mesh_simplify_utils import merge_coplanar_faces, drop_sliver_components

def simplify_for_conversion(tm, target_faces=20000):
    """
    Reduce triangle count before the slow BREP sewing step. The
    rib/bridge/hole solids are extruded prisms -- lots of coplanar
    triangles on flat walls that don't need per-triangle resolution in
    the BREP. Falls back to a no-op (with a hint) if the optional
    `fast-simplification` backend isn't installed.
    """
    if tm is None or len(tm.faces) <= target_faces:
        return tm

    current = len(tm.faces)
    target_reduction = 1.0 - (target_faces / current)
    target_reduction = min(max(target_reduction, 0.0), 0.99)  # clamp to valid (0,1) range

    try:
        simplified = tm.simplify_quadric_decimation(target_reduction)
        print(f"Simplified cut_solid: {current} -> {len(simplified.faces)} faces "
              f"(reduction={target_reduction:.3f})")
        if not simplified.is_watertight:
            print("Warning: simplification broke watertightness, using original mesh instead.")
            return tm
        return simplified
    except Exception as e:
        print(f"Simplification skipped ({e}); "
              f"`pip install fast-simplification` for a speedup here.")
        return tm

def trimesh_to_part_solid(tm, tolerance=0.01):
    """
    Convert a trimesh.Trimesh (expected watertight/volume) into a
    FreeCAD Part solid via Mesh -> makeShapeFromMesh -> sew -> Solid.
    """
    if tm is None or len(tm.vertices) == 0:
        return None

    fc_mesh = trimesh_to_freecad(tm)
    if fc_mesh is None:
        return None

    shape = Part.Shape()
    try:
        shape.makeShapeFromMesh(fc_mesh.Topology, tolerance)
    except Exception as e:
        print(f"trimesh_to_part_solid: makeShapeFromMesh failed: {e}")
        return None

    solid = None
    try:
        solid = Part.makeSolid(shape)
    except Exception:
        try:
            solid = Part.Solid(Part.Shell(shape.Faces))
        except Exception as e:
            print(f"trimesh_to_part_solid: could not build solid, "
                  f"using raw sewn shell instead: {e}")
            solid = shape

    if solid is not None and not solid.isValid():
        try:
            solid.fix(tolerance, tolerance, tolerance)
        except Exception:
            pass

    # An invalid or unclosed cut tool does not make wing_shape.cut() raise
    # -- it either returns a Null shape or sends OCC down its very slow
    # non-solid path. Report it here, where the cause is, instead of
    # letting it surface later as an unexplained failure or a hang.
    if solid is not None:
        print(f"  cut_part: {len(solid.Faces)} faces, {len(solid.Solids)} solid(s), "
              f"valid={solid.isValid()}, closed={solid.isClosed()}")
        if not solid.isValid() or not solid.isClosed():
            print("  WARNING: cut_part is not a valid closed solid -- the wing "
                  "cut will be unreliable and may be extremely slow.")

    return solid

def assemble_final_wing_trimesh(wing_mesh_tm, rib_solid_tm, bridge_solid_tm,
                                 hole_solid_tm, boolean_engine='manifold',
                                 max_fragment_volume=0.0):
    """
    Do the whole assembly in trimesh: wing - (rib - bridge - hole).

    Use this whenever `wing_mesh_tm` is watertight. Both operands are
    already valid volumes by construction (rib/bridge/hole are each
    watertight and manifold booleans preserve that), so there is nothing
    BREP can add here -- and plenty it takes away. The BREP route must
    first sew the cut tool into a Part solid, which loses faces to the
    sewing tolerance (37578 -> 37539, valid=False on this wing) and then
    runs a triangle-soup boolean against a 168k-face wing. Measured on
    this geometry: manifold does the equivalent cut in ~0.3s, OCC takes
    hours. The wing solids here are mesh-derived, so BREP is carrying all
    the cost of exact surface algebra with none of the benefit.

    Returns (final_mesh, cut_solid); final_mesh is None on failure.
    """
    if wing_mesh_tm is None or not wing_mesh_tm.is_watertight:
        print("assemble_final_wing_trimesh: wing mesh is not watertight -- "
              "refusing (the result would have no real internal structure).")
        return None, None
    if rib_solid_tm is None or len(rib_solid_tm.vertices) == 0:
        print("assemble_final_wing_trimesh: no rib_solid, nothing to cut.")
        return None, None

    cut_solid_tm = rib_solid_tm
    for label, other in (("bridges", bridge_solid_tm), ("holes", hole_solid_tm)):
        if other is None or len(other.vertices) == 0:
            continue
        try:
            cut_solid_tm = trimesh.boolean.difference(
                [cut_solid_tm, other], engine=boolean_engine)
            print(f"  rib - {label}: {len(cut_solid_tm.faces)} faces, "
                  f"watertight={cut_solid_tm.is_watertight}")
        except Exception as e:
            print(f"assemble_final_wing_trimesh: rib - {label} failed: {e}")
            return None, cut_solid_tm

    try:
        final = trimesh.boolean.difference(
            [wing_mesh_tm, cut_solid_tm], engine=boolean_engine)
    except Exception as e:
        print(f"assemble_final_wing_trimesh: wing - cut_solid failed: {e}")
        return None, cut_solid_tm

    final = drop_sliver_components(
        final, max_fragment_volume=max_fragment_volume)

    removed = wing_mesh_tm.volume - final.volume
    print(f"Final (trimesh): {len(final.faces)} faces, "
          f"watertight={final.is_watertight}, volume={final.volume:.1f}")
    print(f"  material removed by ribs: {removed:.1f} mm^3 "
          f"({100.0 * removed / wing_mesh_tm.volume:.2f}% of the wing)")
    if removed <= 0:
        print("  WARNING: nothing was removed -- the cut did not take effect.")
    return final, cut_solid_tm


def assemble_final_wing_freecad(wing_shape, rib_solid_tm, bridge_solid_tm, hole_solid_tm,
                                 doc=None, vis=False, mesh_tolerance=0.01,
                                 simplify_target_faces=20000, boolean_engine='manifold'):
    """
    Returns a Part.Shape (the final cut wing) if vis=True and the BREP
    cut succeeds; otherwise returns None but still saves the combined
    cut_solid mesh into the document.
    """
    if rib_solid_tm is None or len(rib_solid_tm.vertices) == 0:
        print("assemble_final_wing_freecad: no rib_solid, nothing to cut.")
        return None

    # ---- Step 1: rib - bridge - hole, entirely in trimesh ----
    cut_solid_tm = rib_solid_tm

    if bridge_solid_tm is not None and len(bridge_solid_tm.vertices) > 0:
        try:
            cut_solid_tm = trimesh.boolean.difference(
                [cut_solid_tm, bridge_solid_tm], engine=boolean_engine)
        except Exception as e:
            print(f"assemble_final_wing_freecad: rib - bridge (trimesh) failed: {e}")
    if cut_solid_tm is None or len(cut_solid_tm.vertices) == 0:
        print("assemble_final_wing_freecad: cut_solid empty after subtracting bridges.")
        return None

    if hole_solid_tm is not None and len(hole_solid_tm.vertices) > 0:
        try:
            cut_solid_tm = trimesh.boolean.difference(
                [cut_solid_tm, hole_solid_tm], engine=boolean_engine)
        except Exception as e:
            print(f"assemble_final_wing_freecad: rib - hole (trimesh) failed: {e}")
    if cut_solid_tm is None or len(cut_solid_tm.vertices) == 0:
        print("assemble_final_wing_freecad: cut_solid empty after subtracting holes.")
        return None

    print(f"cut_solid (trimesh): {len(cut_solid_tm.vertices)} verts, "
          f"{len(cut_solid_tm.faces)} faces, watertight={cut_solid_tm.is_watertight}")

    # ---- Step 2a: coplanar merge -- ONLY if it stays watertight ----
    # cut_solid arrives watertight (rib/bridge/hole are each watertight and
    # manifold booleans preserve that), so the tool is already valid by
    # construction -- nothing here needs repairing. The merge is the only
    # step that can invalidate it: it re-triangulates each planar patch
    # from its merged boundary, so new edges span several original ones
    # and neighbouring faces land mid-edge (T-junctions), leaving
    # unmatched edges. A leaky mesh sews into a shell, Part.makeSolid()
    # then fails ("Creation of solid failed"), and cutting the wing with
    # the raw shell that gets used instead either returns a Null shape or
    # grinds for hours in OCC's non-solid path.
    #
    # Measured on WingR2 (same geometry, merge the only difference):
    #   without merge: 1560 faces, watertight -> tool valid=True,
    #                  closed=True  -> cut 13.0s, result valid=True
    #   with merge   : 1431 faces, leaky      -> tool valid=False
    #                  -> cut 12.5s, result valid=False
    # Same runtime and the same 509.3 mm^3 removed, so the merge's only
    # measurable effect here is to destroy validity. Keep it when it is
    # genuinely lossless, discard it otherwise.
    pre_merge_tm = cut_solid_tm
    merged_tm = merge_coplanar_faces(cut_solid_tm)
    if merged_tm is not None and merged_tm.is_watertight:
        cut_solid_tm = merged_tm
    else:
        print(f"Coplanar merge broke watertightness -- discarding it and "
              f"keeping the pre-merge mesh ({len(pre_merge_tm.faces)} faces), "
              f"which is still a valid volume.")
        cut_solid_tm = pre_merge_tm

    # ---- Step 2b: lossy simplification -- SKIP unless coplanar merge
    # wasn't enough. Quadric decimation collapses small notches (like the
    # bridge_height-scale cavities cut_solid has wherever a bridge/hole
    # was subtracted) just as readily as real redundancy, since their
    # geometric error is small relative to the whole mesh. That's what
    # was silently erasing your bridges/holes from the exported STL. ----
    print(f"cut_solid after coplanar merge: {len(cut_solid_tm.faces)} faces")
    if len(cut_solid_tm.faces) > simplify_target_faces:
        print("Still large after coplanar merge -- applying quadric decimation "
              "(this CAN remove thin features; verify the exported STL carefully).")
        cut_solid_tm = simplify_for_conversion(cut_solid_tm, simplify_target_faces)
    else:
        print("Skipping quadric decimation -- coplanar merge alone was sufficient.")

    # ---- Add the (simplified) cut_solid mesh to the project and save,
    # regardless of whether we go on to do the full wing cut. This is
    # the useful, fast-to-get-to result -- worth persisting even if the
    # slow BREP step below is skipped or fails. ----
    if doc is not None:
        try:
            from viz_utils import show_mesh, fit_view
            fc_mesh = trimesh_to_freecad(cut_solid_tm)
            if fc_mesh:
                show_mesh(fc_mesh, doc, "CutSolid", color=(0.9, 0.4, 0.1), transparency=20)
                doc.recompute()
                print("Cut solid (rib - bridge - hole) visualised.")
        except Exception as e:
            print(f"Visualisation error: {e}")
        doc.save()
        print("Document saved (cut_solid mesh).")
        fit_view(doc)

    # ---- Only do the expensive wing - cut_solid BREP boolean if requested ----
    if not vis:
        print("vis_final is False -- skipping full wing cut; cut_solid mesh saved above.")
        return None

    # ---- Step 3: convert ONLY this one mesh to a Part solid ----
    print("Converting combined cut_solid to Part solid (single conversion)...")
    cut_part = trimesh_to_part_solid(cut_solid_tm, mesh_tolerance)
    if cut_part is None:
        print("assemble_final_wing_freecad: cut_solid conversion failed, aborting.")
        return None

    # ---- Step 4: the one real BREP boolean, against the actual wing ----
    print("Performing final wing - cut_solid boolean (BREP, on wing_shape)...")
    try:
        wing_obj = wing_shape.cut(cut_part)
    except Exception as e:
        print(f"assemble_final_wing_freecad: wing - cut_solid failed: {e}")
        return None

    print(f"Final wing_obj: {len(wing_obj.Faces)} faces, "
          f"{len(wing_obj.Solids)} solid(s), valid={wing_obj.isValid()}")

    try:
        from viz_utils import show_shape, fit_view
        show_shape("WingFinal", wing_obj, doc, color=(0.75, 0.75, 0.75), transparency=10)
        doc.recompute()
        doc.save()
        print("Final wing_obj visualised and saved.")
        fit_view(doc)
    except Exception as e:
        print(f"Visualisation error: {e}")

    return wing_obj
