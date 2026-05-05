# point_utils.py
import FreeCAD
import Part

def collect_midpoints_from_wing_and_ribs(wing_shape, rib_solids, z_min, z_max, z_step):
    """
    Slice the original wing horizontally. For each rib solid,
    intersect it with the same plane, take the midpoint of the intersection edge.
    Returns dict {rib_index: list_of_points}
    """
    points_by_rib = {i: [] for i in range(len(rib_solids))}
    z = z_min
    slice_count = 0
    while z <= z_max:
        slice_count += 1
        # Small epsilon to increase chance of hitting thin ribs
        eps = 0.01
        for idx, rib in enumerate(rib_solids):
            # Try both exact z and z+eps (and z-eps for safety)
            for zz in (z, z + eps, z - eps):
                try:
                    slice_compound = rib.slice(FreeCAD.Vector(0,0,1), zz)
                    if slice_compound and hasattr(slice_compound, 'Edges') and slice_compound.Edges:
                        edge = slice_compound.Edges[0]
                        if len(edge.Vertexes) >= 2:
                            p1 = edge.Vertexes[0].Point
                            p2 = edge.Vertexes[-1].Point
                            mid = (p1 + p2) * 0.5
                            points_by_rib[idx].append(mid)
                            break   # stop after first successful intersection at this z
                except Exception:
                    continue
        z += z_step

    total = sum(len(pts) for pts in points_by_rib.values())
    print(f"Collected {total} midpoints over {slice_count} slices.")
    return points_by_rib

def show_points_per_rib(points_by_rib, doc, prefix="RibPoints"):
    """Creates a separate Part::Feature with vertices for each rib."""
    count = 0
    for idx, pts in points_by_rib.items():
        if not pts:
            continue
        vertices = [Part.Vertex(p) for p in pts]
        compound = Part.Compound(vertices)
        obj = doc.addObject("Part::Feature", f"{prefix}_{idx}")
        obj.Shape = compound
        count += 1
    doc.recompute()
    total = sum(len(p) for p in points_by_rib.values())
    print(f"Visualized {count} ribs with points (total {total} points).")