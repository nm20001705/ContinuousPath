import FreeCAD
import Part
import Free
from PySide import QtGui
from PySide import QtCore

# shape = Part.Shape()
# shape.read(r"C:\Users\natha\git\ContinuousPath\.in\test_wing.step")

# hol=shape.makeThickness([shape.Faces[0]],-0.4,1000)
# Part.show(hol)
s = FreeCADGui.Selection.getSelection()
try:
    shape1=s[0].Shape
except:
    print "Wrong selection"

myObject = App.ActiveDocument.addObject("Part::Feature","Shell")

Shell = shape1.makeThickness([], 3.0, 0.1)
myObject.Shape = Shell