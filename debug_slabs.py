# debug_slabs.py – visualise slab surfaces

import FreeCAD
import FreeCADGui
import Part
import Mesh
import MeshPart
import numpy as np
import trimesh
import tempfile
import os
from slab_utils import create_angled_grid_lines, get_slab_surface, add_slab_to_doc
from viz_utils import fit_view

# --- Load wing (adjust path) ---
doc = FreeCAD.open(r"C:\Users\natha\Desktop\plane\3D\Slop3r V-tail slope glider 1.2 m span - 4647489\0_make_struct\fin.FCStd")
wing_obj = doc.getObject('Extrude001')
wing_shape = wing_obj.Shape

# Convert to trimesh (use existing function from earlier)
def shape_to_trimesh(shape, deflection=0.5, angular=0.5):
    with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp:
        tmp_stl = tmp.name
    try:
        mesh_obj = MeshPart.meshFromShape(shape, LinearDeflection=deflection, AngularDeflection=angular, Relative=False)
        mesh_obj.write(tmp_stl)
        tm = trimesh.load(tmp_stl, force='mesh')
        os.unlink(tmp_stl)
        return tm
    except:
        return None

wing_mesh = shape_to_trimesh(wing_shape)
if wing_mesh is None:
    raise RuntimeError("Failed to convert wing to mesh")

bb = wing_shape.BoundBox

# --- Create parameters object (minimal) ---
class DummyParams:
    pass
params = DummyParams()
params.plane_normal = np.array([0,1,0])   # XZ plane
params.plane_axis_u = np.array([1,0,0])
params.plane_axis_v = np.array([0,0,1])
params.rib_spacing = 20.0
params.rib_angle = 30.0
params.grid_orientation = 0.0
params.primary_dir = np.array([0,0,1])  # or None
params.rib_width = 0.2   # not used for surfaces

# --- Generate grid lines ---
lines1, lines2 = create_angled_grid_lines(bb, params)
all_lines = lines1 + lines2
print(f"Generated {len(all_lines)} rib lines")

# --- For each line, create slab surface and display ---
slab_objects = []
for idx, (start, end) in enumerate(all_lines):
    face = get_slab_surface(wing_mesh, start, end, params.plane_normal)
    if face is not None:
        obj = add_slab_to_doc(face, doc, f"Slab_{idx}", color=(0.5,0.7,0.9), transparency=30)
        slab_objects.append(obj)
        print(f"Added slab {idx}")
    else:
        print(f"Slab {idx} failed")

doc.recompute()
fit_view(doc)
print("Done. Check the FreeCAD view.")