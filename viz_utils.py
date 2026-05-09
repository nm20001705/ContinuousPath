# viz_utils.py
import FreeCAD
import FreeCADGui
import Part

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