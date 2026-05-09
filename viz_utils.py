# viz_utils.py
import FreeCAD
import FreeCADGui
import Part
import math

def show_shape(obj_name, shape, doc, color=None, transparency=None,
               line_color=None, line_width=None, point_size=None,
               display_mode=None):
    """
    Create a Part::Feature object and apply view properties if GUI is active.
    """
    if not shape or shape.isNull():
        print(f"Warning: Cannot show {obj_name} – shape is null.")
        return None
    obj = doc.addObject("Part::Feature", obj_name)
    obj.Shape = shape
    if FreeCAD.GuiUp:
        vp = FreeCADGui.ActiveDocument.getObject(obj.Name)
        if vp:
            if color is not None:
                vp.ViewObject.ShapeColor = color
            if transparency is not None:
                vp.ViewObject.Transparency = min(max(transparency, 0), 100)  # clamp 0-100
            if line_color is not None:
                vp.ViewObject.LineColor = line_color
            if line_width is not None:
                vp.ViewObject.LineWidth = line_width
            if point_size is not None:
                vp.ViewObject.PointSize = point_size
            if display_mode is not None:
                vp.ViewObject.DisplayMode = display_mode
            else:
                # Ensure default display mode is "Shaded" to show transparency
                if transparency is not None:
                    vp.ViewObject.DisplayMode = "Shaded"
    doc.recompute()
    return obj

def show_cut_wing(shape, doc, transparency=80):
    return show_shape("CutWing", shape, doc, transparency=transparency, display_mode="Shaded")

def show_rib_centre_lines(lines, doc, line_color=(1.0,0.0,0.0)):
    compound = Part.Compound([line for line in lines])
    return show_shape("RibCentreLines", compound, doc, line_color=line_color, line_width=2)

def show_midpoints(points, doc, point_size=5, color=(1.0,0.0,0.0)):
    vertices = [Part.Vertex(p) for p in points]
    compound = Part.Compound(vertices)
    return show_shape("Midpoints", compound, doc, color=color, point_size=point_size, display_mode="Points")

def show_rib_wires(wires, doc, line_color=(0.0,1.0,0.0), line_width=2):
    compound = Part.Compound(wires)
    return show_shape("RibWires", compound, doc, line_color=line_color, line_width=line_width)

def show_bridges(shape, doc, color=(0.0,0.8,0.0), transparency=30):
    return show_shape("Bridges", shape, doc, color=color, transparency=transparency, display_mode="Shaded")

def show_final_solid(shape, doc, color=(0.8,0.8,0.8), transparency=80):
    return show_shape("WingWithBridges", shape, doc, color=color, transparency=transparency, display_mode="Shaded")

def fit_view(doc):
    if FreeCAD.GuiUp:
        FreeCADGui.SendMsgToActiveView("ViewFit")

def show_rib_centre_surfaces(faces, doc, color=(0.6, 0.6, 0.9), transparency=60):
    """
    Display the rib centre surfaces (mid‑plane trimmed to wing).
    """
    if not faces:
        return None
    compound = Part.Compound(faces)
    return show_shape("RibCentreSurfaces", compound, doc, color=color, transparency=transparency)

def show_rib_solids_trimmed(rib_solids, doc, color=(0.5, 0.7, 0.5), transparency=50):
    """
    Display the rib solids that have been trimmed to the wing.
    """
    if not rib_solids:
        return None
    compound = Part.Compound(rib_solids)
    return show_shape("RibSolids", compound, doc, color=color, transparency=transparency)

def show_rib_centre_edges(edges, doc, line_color=(0.2, 0.5, 1.0), line_width=2):
    """
    Display the intersection edges of the rib planes with the wing.
    """
    if not edges:
        return None
    compound = Part.Compound(edges)
    return show_shape("RibCentreEdges", compound, doc, line_color=line_color, line_width=line_width)

def show_rib_segments(pieces, doc, color=(0.2, 0.8, 0.4), transparency=30):
    """
    Display all rib face segments as a single compound.
    pieces: flat list of Part.Face from split_rib_faces_by_crossings.
    """
    if not pieces:
        print("show_rib_segments: no pieces to show.")
        return None
    valid = [p for p in pieces if p is not None and not p.isNull()]
    if not valid:
        print("show_rib_segments: all pieces are null.")
        return None
    compound = Part.Compound(valid)
    return show_shape("RibSegments", compound, doc,
                      color=color, transparency=transparency)

def add_ellipse_to_doc(center_3d, a, b, ang, normal, doc, color=(1.0,0.0,0.0)):
    """Create an ellipse wire and face in 3D."""
    # We need to build the ellipse in the local plane and rotate it.
    # For simplicity, create a circle and scale if a != b.
    # But proper ellipse: define a wire of 3D points.
    points = []
    for t in range(0, 360, 10):
        rad = math.radians(t)
        # Local coordinates
        x = a * math.cos(rad)
        y = b * math.sin(rad)
        # Rotate by ang around normal
        # Build rotation matrix: X' = (cos ang, sin ang), Y' = (-sin ang, cos ang) in the plane
        # Need two axes in the plane perpendicular to normal.
        # We assume centre_3d is in 3D; we need to define u and v orthonormal in the plane.
        # Use arbitrary u and v from face's normal.
        n = normal.normalize()
        if abs(n.x) < 0.9:
            u = FreeCAD.Vector(1,0,0).cross(n).normalize()
        else:
            u = FreeCAD.Vector(0,1,0).cross(n).normalize()
        v = n.cross(u).normalize()
        # Rotate (x,y) by angle ang
        xx = x * math.cos(ang) - y * math.sin(ang)
        yy = x * math.sin(ang) + y * math.cos(ang)
        pt_3d = center_3d + u * xx + v * yy
        points.append(pt_3d)
    wire = Part.makePolygon(points + [points[0]])
    face = Part.Face(wire)
    obj = doc.addObject("Part::Feature", "Ellipse")
    obj.Shape = face
    if FreeCAD.GuiUp:
        obj.ViewObject.ShapeColor = color
        obj.ViewObject.Transparency = 30
    doc.recompute()
    return obj
