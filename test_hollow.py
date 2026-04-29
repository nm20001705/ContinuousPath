import pyvista as pv
import numpy as np

def hollow_mesh(path, thickness=0.4):

    mesh = pv.read(path)

    # clean surface
    mesh = mesh.extract_surface().triangulate().clean()

    # 1. create sampling grid (IMPORTANT)
    bounds = mesh.bounds
    spacing = thickness / 3.0

    grid = pv.UniformGrid()

    nx = int((bounds[1] - bounds[0]) / spacing)
    ny = int((bounds[3] - bounds[2]) / spacing)
    nz = int((bounds[5] - bounds[4]) / spacing)

    grid.dimensions = (nx, ny, nz)

    grid.origin = (bounds[0], bounds[2], bounds[4])
    grid.spacing = (spacing, spacing, spacing)

    # 2. compute signed distance to surface
    sdf = grid.compute_implicit_distance(mesh)

    # 3. extract inner surface (offset ≈ thickness)
    inner = sdf.contour([thickness])

    # 4. boolean difference
    shell = mesh.boolean_difference(inner)

    return shell


shell = hollow_mesh(".in/test_wing.stl", thickness=0.4)
shell.save(".out/test_hollow.stl")