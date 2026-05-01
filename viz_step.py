import cadquery as cq

shape = cq.importers.importStep(".in/hollow.step")
solid = shape.solids().val()