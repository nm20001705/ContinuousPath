import numpy as np
import trimesh
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import trimesh
from scipy.ndimage import binary_erosion

class OuterGeometryProvider:
    def get_contour(self, z: float) -> np.ndarray:
        raise NotImplementedError

class CustomShape(OuterGeometryProvider):
    def __init__(self, points):
        self.points = np.array(points)

        # assume y-axis is "height"
        self.y_min = np.min(self.points[:, 1])
        self.y_max = np.max(self.points[:, 1])
        self.height = self.y_max - self.y_min

    def get_contour(self, z):
        return self.points

    def get_contour(self, z):
        return self.points

class GridBuilder(OuterGeometryProvider):
    def __init__(self):
        pass

class Path:
    def __init__(self, path):
        self.path = np.array(path)

        # ordered unique heights (critical fix)
        self.heights = np.unique(self.path[:, -1])
        self.n_layers = len(self.heights)

    def get_layer(self, n):
        height = self.heights[n]
        layer = self.path[self.path[:, -1] == height]
        return Path(layer)
    
# =========================================================
# PARAMETER CONTAINER
# =========================================================
class ContourParams:
    def __init__(self,
                 layer_height=0.3,
                 line_width=0.8,
                 rib_clearance=0.5,
                 grid_angle=45,
                 grid_spacing=5,
                 grid_orientation="x"):  # "x", "y", "both"

        self.layer_height = layer_height
        self.line_width = line_width
        self.rib_clearance = rib_clearance
        self.grid_angle = np.deg2rad(grid_angle)
        self.grid_spacing = grid_spacing
        self.grid_orientation = grid_orientation
        self.radius = line_width/2

# =========================================================
# GRID GENERATOR
# =========================================================
class GridGenerator:
    def __init__(self, params: ContourParams):
        self.p = params

    def generate_ribs(self, outer):
        ribs = []

        bbox_min = outer.min(axis=0)
        bbox_max = outer.max(axis=0)

        spacing = self.p.grid_spacing
        angle = self.p.grid_angle

        # generate parallel lines
        direction = np.array([np.cos(angle), np.sin(angle)])
        normal = np.array([-direction[1], direction[0]])

        extent = np.linalg.norm(bbox_max - bbox_min)

        offsets = np.arange(-extent, extent, spacing)

        for o in offsets:
            origin = bbox_min + o * normal

            # long line across domain
            p0 = origin - 2 * extent * direction
            p1 = origin + 2 * extent * direction

            ribs.append((p0, p1))

        return ribs

# =========================================================
# PATH BUILDER
# =========================================================
class PathBuilder:
    def __init__(self, params: ContourParams):
        self.p = params

    def build_layer(self, outer, ribs):
        path = []

        for i in range(len(outer) - 1):
            p0 = outer[i]
            p1 = outer[i + 1]

            path.append(p0)

            for (r0, r1) in ribs:
                # simple intersection heuristic
                if self._near_segment(r0, p0, p1):
                    excursion = self._make_excursion(r0, r1, outer)
                    path.extend(excursion)

        path.append(outer[-1])
        return np.array(path)

    def _near_segment(self, pt, a, b, tol=1.0):
        # distance from point to segment
        ap = pt - a
        ab = b - a
        t = np.dot(ap, ab) / np.dot(ab, ab)
        t = np.clip(t, 0, 1)
        closest = a + t * ab
        return np.linalg.norm(pt - closest) < tol

    def _make_excursion(self, r0, r1, outer):
        # shorten rib so it doesn't hit outer wall
        center = (r0 + r1) / 2
        direction = (r1 - r0)
        direction /= np.linalg.norm(direction)

        length = 5.0  # simple fixed length for now

        p_in = center - direction * length
        p_out = center + direction * length

        return [p_in, p_out, p_in]

# =========================================================
# WING GENERATOR
# =========================================================
class ToolpathGenerator:
    def __init__(self, params, geometry: OuterGeometryProvider):
        self.p = params
        self.geometry = geometry
        self.grid = GridGenerator(params)
        self.path_builder = PathBuilder(params)

    def build_path(self):
        layers = int(self.geometry.height / self.p.layer_height)
        all_points = []

        for i in range(layers):
            z = i * self.p.layer_height

            outer = self.geometry.get_contour(z)

            ribs = self.grid.generate_ribs(outer)
            layer_path = self.path_builder.build_layer(outer, ribs)

            layer_3d = np.column_stack([
                layer_path[:, 0],
                layer_path[:, 1],
                np.full(len(layer_path), z)
            ])

            all_points.append(layer_3d)

        return Path(np.vstack(all_points))

# =========================================================
# MESH GENERATOR
# =========================================================
class MeshBuilder:
    def __init__(self, radius):
        self.radius = radius

    def build(self, path):
        path = path.path
        tangents = np.gradient(path, axis=0)
        tangents /= np.linalg.norm(tangents, axis=1, keepdims=True)

        ref = np.array([0, 0, 1])
        normals = []
        binormals = []

        for t in tangents:
            n = np.cross(t, ref)
            if np.linalg.norm(n) < 1e-6:
                ref = np.array([0, 1, 0])
                n = np.cross(t, ref)

            n /= np.linalg.norm(n)
            b = np.cross(t, n)

            normals.append(n)
            binormals.append(b)

        normals = np.array(normals)
        binormals = np.array(binormals)

        sections = 10
        circle = np.linspace(0, 2*np.pi, sections, endpoint=False)

        vertices = []
        faces = []

        for i, p in enumerate(path):
            ring = []
            for theta in circle:
                offset = np.cos(theta)*normals[i] + np.sin(theta)*binormals[i]
                ring.append(p + self.radius * offset)
            vertices.append(ring)

        vertices = np.array(vertices)

        for i in range(len(path)-1):
            for j in range(sections):
                a = i*sections + j
                b = i*sections + (j+1)%sections
                c = (i+1)*sections + j
                d = (i+1)*sections + (j+1)%sections

                faces.append([a, c, b])
                faces.append([b, c, d])

        return trimesh.Trimesh(
            vertices=vertices.reshape(-1, 3),
            faces=np.array(faces),
            process=True
        )

class TrimeshGeometry(OuterGeometryProvider):
    def __init__(self, mesh: trimesh.Trimesh, params):
        self.params = params
        self.mesh = mesh
        self.z_min = mesh.bounds[0][2]
        self.z_max = mesh.bounds[1][2]
        self.height = self.z_max - self.z_min

    def get_contour(self, z):
        section = self.mesh.section(
            plane_origin=[0, 0, z],
            plane_normal=[0, 0, 1]
        )

        if section is None:
            return np.empty((0, 2))

        # convert to 2D planar coordinates
        slice_2D, _ = section.to_planar()

        # extract discrete paths
        paths = slice_2D.discrete

        if len(paths) == 0:
            return np.empty((0, 2))

        # pick largest loop (outer boundary)
        largest = max(paths, key=lambda p: self._polygon_area(p))

        return np.array(largest)

    def _polygon_area(self, pts):
        x = pts[:, 0]
        y = pts[:, 1]
        return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

    def hollow_mesh(self):
        pitch = self.params.line_width

        mesh = self.mesh.copy()
        mesh = mesh.subdivide_to_size(max_edge=pitch * 2)

        vox = mesh.voxelized(pitch)
        filled = vox.fill()

        # --- KEY FIX HERE ---
        matrix = filled.matrix

        # erode inward (controls wall thickness)
        eroded = binary_erosion(matrix, iterations=1)

        shell = matrix ^ eroded  # XOR → shell

        # rebuild voxel grid
        shell_vox = trimesh.voxel.VoxelGrid(shell, transform=filled.transform)

        shell_mesh = shell_vox.marching_cubes

        return TrimeshGeometry(shell_mesh, self.params)

# =========================================================
# MAIN
# =========================================================

def plot_path(path, start=0, step=1, stop=None, duration=3):
    if stop is None:
        stop = path.n_layers

    layers = list(range(start, stop, step))

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    line, = ax.plot([], [], lw=1)

    # set stable bounds (important)
    all_pts = path.path
    xmin = all_pts[:, 0].min()
    xmax = all_pts[:, 0].max()
    width = xmax-xmin
    ymin = all_pts[:, 1].min()
    ymax = all_pts[:, 1].max()
    height = ymax-ymin
    ax.set_xlim(all_pts[:, 0].min()-0.9*width, all_pts[:, 0].max()+1.1*width)
    ax.set_ylim(all_pts[:, 1].min()-0.9*height, all_pts[:, 1].max()+1.1*height)

    def update(i):
        layer = path.get_layer(layers[i]).path

        x = layer[:, 0]
        y = layer[:, 1]

        line.set_data(x, y)
        return line,

    fps = 30
    frames = len(layers)
    interval = 1000 * duration / frames

    ani = FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=interval,
        blit=True
    )

    plt.show()
    return ani

if __name__ == "__main__":
    params = ContourParams(layer_height=0.1, line_width=0.4, rib_clearance=0.0, grid_angle=45, grid_spacing=5, grid_orientation='z')
    mesh_stl = trimesh.load(".in/test_wing.stl")
    geometry = TrimeshGeometry(mesh_stl, params=params).hollow_mesh()

    geometry.mesh.export(".out/test.stl")

    pass