import FreeCAD
import Part
import math
import numpy as np

# ============================================================
# PARAMETERS
# ============================================================
class LWInfillParams:
    def __init__(self,
                 nozzle_diameter    = 0.4,
                 wall_thickness     = None,
                 rib_spacing        = 5.0,
                 rib_width          = None,
                 split_angle        = 45.0,
                 chord_curves       = None,
                 guide_rails        = None,
                 rib_angle          = 30.0):   # degrees, half-angle between rib families
        self.nozzle_diameter = nozzle_diameter
        self.wall_thickness = wall_thickness or nozzle_diameter
        self.rib_spacing = rib_spacing
        self.rib_width = rib_width or (2.0 * nozzle_diameter)
        self.split_angle = split_angle
        self.chord_curves = chord_curves
        self.guide_rails = guide_rails
        self.rib_angle = rib_angle


# ============================================================
# FACE SELECTION (automatic)
# ============================================================
def find_bottom_face(body):
    """
    Return the face whose centre of mass has the smallest Z coordinate.
    (Usually the flat bottom of the wing, best for drawing the grid.)
    """
    best_face = None
    best_z = float('inf')
    for face in body.Faces:
        try:
            com = face.CenterOfMass
            if com.z < best_z:
                best_z = com.z
                best_face = face
        except Exception:
            pass
    return best_face


# ============================================================
# GRID LINE GENERATION
# ============================================================
def create_angled_grid_lines(face, spacing, angle_deg):
    """
    Generate two families of parallel lines on the face,
    at +angle_deg and -angle_deg relative to the face's long axis.
    Returns (lines_family1, lines_family2).
    """
    bb = face.BoundBox
    z = bb.ZMin   # work on the horizontal plane of the face

    # Determine the primary axis (longer side of the bounding box)
    if bb.XLength > bb.YLength:
        primary_dir = FreeCAD.Vector(1, 0, 0)   # X is span
    else:
        primary_dir = FreeCAD.Vector(0, 1, 0)   # Y is span

    # Perpendicular horizontal vector
    perp_primary = FreeCAD.Vector(-primary_dir.y, primary_dir.x, 0)

    # Direction vectors for the two rib families
    angle_rad = math.radians(angle_deg)
    v1 = primary_dir * math.cos(angle_rad) + perp_primary * math.sin(angle_rad)
    v2 = primary_dir * math.cos(angle_rad) - perp_primary * math.sin(angle_rad)

    # Bounding box corners in 2D (XY plane at z)
    corners_2d = [
        FreeCAD.Vector(bb.XMin, bb.YMin, z),
        FreeCAD.Vector(bb.XMax, bb.YMin, z),
        FreeCAD.Vector(bb.XMin, bb.YMax, z),
        FreeCAD.Vector(bb.XMax, bb.YMax, z)
    ]

    def generate_lines(v):
        # Perpendicular direction in XY plane
        perp = FreeCAD.Vector(-v.y, v.x, 0)
        # Project all corners onto perp to get the range
        proj_vals = [c.dot(perp) for c in corners_2d]
        min_proj = min(proj_vals)
        max_proj = max(proj_vals)

        # Number of lines needed
        num_lines = int((max_proj - min_proj) / spacing) + 3
        center_pt = FreeCAD.Vector(bb.Center.x, bb.Center.y, z)
        line_len = max(bb.XLength, bb.YLength) * 3   # generous over‑coverage

        lines = []
        for i in range(-1, num_lines + 1):
            offset = min_proj + i * spacing
            # Move the center point perpendicularly to the correct offset
            p0 = center_pt + perp * (offset - center_pt.dot(perp))
            line_start = p0 - v * line_len
            line_end   = p0 + v * line_len
            lines.append(Part.makeLine(line_start, line_end))
        return lines

    return generate_lines(v1), generate_lines(v2)


# ============================================================
# CREATE THIN RIB FACES FROM LINES
# ============================================================
def create_rib_faces(lines, plane_normal, rib_width):
    """
    Convert each line into a thin rectangular face oriented perpendicular
    to the plane (i.e., vertical for vase printing).
    """
    half_w = rib_width / 2.0
    faces = []
    for line in lines:
        try:
            start = line.Vertexes[0].Point
            end = line.Vertexes[-1].Point
            dir_vec = (end - start).normalize()
            # A vector perpendicular to the line AND to the plane normal (so the face is vertical)
            perp = dir_vec.cross(plane_normal).normalize() * half_w
            p1 = start + perp
            p2 = start - perp
            p3 = end - perp
            p4 = end + perp
            wire = Part.makePolygon([p1, p2, p3, p4, p1])
            faces.append(Part.Face(wire))
        except Exception:
            pass
    return faces


# ============================================================
# MAIN GENERATION (without split)
# ============================================================
def generate_lw_infill(body, params, doc=None):
    if doc is None:
        doc = FreeCAD.ActiveDocument

    # 1. Face selection
    work_face = find_bottom_face(body)
    if not work_face:
        raise ValueError("No valid face found on the body.")
    face_normal = work_face.normalAt(0.5, 0.5)

    # 2. Create two families of grid lines
    print(f"Generating grid lines on face (area={work_face.Area:.1f} mm²)...")
    lines1, lines2 = create_angled_grid_lines(work_face, params.rib_spacing, params.rib_angle)
    print(f"  Family1 lines: {len(lines1)}, Family2 lines: {len(lines2)}")

    # 3. Convert lines to rib faces
    rib_faces1 = create_rib_faces(lines1, face_normal, params.rib_width)
    rib_faces2 = create_rib_faces(lines2, face_normal, params.rib_width)
    print(f"  Rib faces: {len(rib_faces1)} + {len(rib_faces2)}")

    if not rib_faces1 and not rib_faces2:
        raise ValueError("No rib faces generated.")

    bbox = body.BoundBox
    z_extent = bbox.ZMax - bbox.ZMin
    extrude_dist = z_extent * 2.0
    extrude_dir = FreeCAD.Vector(0, 0, 1)

    # Helper: extrude a list of faces
    def extrude_faces(faces):
        solids = []
        for face in faces:
            try:
                rib = face.extrude(extrude_dir * extrude_dist)
                rib.translate(FreeCAD.Vector(0, 0, bbox.ZMin - z_extent * 0.5))
                solids.append(rib)
            except Exception:
                pass
        return solids

    solids1 = extrude_faces(rib_faces1)
    solids2 = extrude_faces(rib_faces2)

    # Helper: fuse a list of solids
    def fuse_solids(solids):
        if not solids:
            return None
        fused = solids[0]
        for s in solids[1:]:
            try:
                fused = fused.fuse(s)
            except Exception:
                pass
        return fused

    fused1 = fuse_solids(solids1)
    fused2 = fuse_solids(solids2)

    # Combine the two families (no extra tilt needed – the ribs are already at ±rib_angle)
    rib_parts = [p for p in [fused1, fused2] if p is not None]
    if not rib_parts:
        raise ValueError("No rib solids to cut with.")
    fused_ribs = rib_parts[0]
    for p in rib_parts[1:]:
        try:
            fused_ribs = fused_ribs.fuse(p)
        except Exception:
            pass

    print("Cutting internal structure...")
    result = body
    try:
        result = body.cut(fused_ribs)
    except Exception as e:
        print(f"Cut failed: {e}")
    return result

    def tilt_ribs(ribs, angle):
        if angle == 0.0 or ribs is None:
            return ribs
        rot_mat = FreeCAD.Matrix()
        rot_mat.rotateY(angle)
        T_to_origin = FreeCAD.Matrix()
        T_to_origin.move(-rot_center)
        T_back = FreeCAD.Matrix()
        T_back.move(rot_center)
        final_mat = T_back * rot_mat * T_to_origin
        return ribs.transformGeometry(final_mat)

    if tilt != 0.0:
        print(f"Tilting ribs: ±{params.tilt_angle:.1f}°")
    fused1 = tilt_ribs(fused1, tilt)   # +angle
    fused2 = tilt_ribs(fused2, -tilt)  # -angle

    # Combine both families
    rib_parts = [p for p in [fused1, fused2] if p is not None]
    if not rib_parts:
        raise ValueError("No rib solids to cut with.")
    fused_ribs = rib_parts[0]
    for p in rib_parts[1:]:
        try:
            fused_ribs = fused_ribs.fuse(p)
        except Exception:
            pass

    # Cut
    print("Cutting internal structure...")
    result = body
    try:
        result = body.cut(fused_ribs)
    except Exception as e:
        print(f"Cut failed: {e}")
    return result

# ============================================================
# EXAMPLE USAGE
# ============================================================
if __name__ == "__main__":
    # Open your wing file
    doc_path = r"C:\Users\natha\git\ContinuousPath\wing.FCStd"
    doc = FreeCAD.open(doc_path)

    # Get the wing body (assuming the Pad object is the solid)
    wing = doc.getObject("Pad")
    if not wing:
        raise RuntimeError("Object 'Pad' not found in document.")

    # Set parameters
    params = LWInfillParams(
        nozzle_diameter = 0.4,
        rib_spacing     = 10.0,
        rib_width       = 0.1,
        rib_angle       = 30.0   # +30° and -30° from span
    )

    # Generate the infill (no chord curves yet)
    result_body = generate_lw_infill(wing.Shape, params, doc)

    # Add result to document
    obj = doc.addObject("Part::Feature", "WingWithInfill")
    obj.Shape = result_body
    doc.recompute()
    print("Lightweight wing created successfully.")
    