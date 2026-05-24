import pyvista as pv
import numpy as np
from copy import deepcopy
import tqdm
import warnings
from itertools import permutations
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
from svv.utils.remeshing.remesh import remesh_surface_2d, remesh_volume
from svv.utils.pyvista_compat import compute_cell_quality
import tetgen


class BoundaryLayer(object):
    def __init__(self, surface, layer_thickness=1.0, layer_thickness_ratio=0.5,
                 constant_thickness=True, max_layer_thickness=None, min_layer_thickness=None,
                 number_of_sublayers=3, number_of_substeps=100, sublayer_ratio=0.3,
                 relaxation=0.01, local_correction_factor=0.45, include_surface_cells=True,
                 include_sidewall_cells=True, negate_warp_vectors=True, volume_cell_entity_id=0,
                 inner_surface_cell_entity_id=1, outer_surface_cell_entity_id=2,
                 use_warp_vector_magnitude_as_thickness=False, layer_thickness_array=None,
                 warp_vector_array="Normals", remesh_vol=False, combine=True):
        self.surface = surface
        self.warp_vector_array = warp_vector_array
        self.surface_cells = {}
        self.surface_cell_array = self.surface.faces.reshape(-1, 4)[:, 1:]
        self.point_cells = {}
        for i in range(surface.n_points):
            self.point_cells[i] = []
        for i in tqdm.trange(surface.n_cells, desc="Processing Surface Cells", leave=False):
            cell = surface.get_cell(i)
            data = {}
            data['points'] = cell.points
            data['point_ids'] = cell.point_ids
            data['type'] = cell.type
            for j in cell.point_ids:
                self.point_cells[j].append(i)
            self.surface_cells[i] = data
        self.point_neighbors = {}
        for i in tqdm.trange(surface.n_points, desc="Processing Point Neighbors", leave=False):
            self.point_neighbors[i] = []
            for j in self.point_cells[i]:
                cell = self.surface_cells[j]
                for k in cell['point_ids']:
                    if k != i:
                        self.point_neighbors[i].append(k)
        if warp_vector_array not in surface.point_data.keys():
            if warp_vector_array == "Normals":
                print("Normals not found, computing normals.")
                self.surface.compute_normals(inplace=True)
            else:
                raise ValueError("Warp vector array not found.")
        self.layer_thickness = layer_thickness
        self.layer_thickness_ratio = layer_thickness_ratio
        self.constant_thickness = constant_thickness
        self.max_layer_thickness = max_layer_thickness
        self.min_layer_thickness = min_layer_thickness
        self.number_of_sublayers = number_of_sublayers
        self.number_of_substeps = number_of_substeps
        self.sublayer_ratio = sublayer_ratio
        self.relaxation = relaxation
        self.local_correction_factor = local_correction_factor
        self.include_surface_cells = include_surface_cells
        self.include_sidewall_cells = include_sidewall_cells
        self.negate_warp_vectors = negate_warp_vectors
        self.boundary_layers = pv.UnstructuredGrid()
        self.volume_cell_entity_id = volume_cell_entity_id
        self.inner_surface_cell_entity_id = inner_surface_cell_entity_id
        self.outer_surface_cell_entity_id = outer_surface_cell_entity_id
        self.use_warp_vector_magnitude_as_thickness = use_warp_vector_magnitude_as_thickness
        self.layer_thickness_array = layer_thickness_array
        self.warp_vectors = None
        self.remesh_vol = remesh_vol
        self.tolerance = 1E-14
        self.layers = None
        self.combine = combine

    def generate(self):
        grid = deepcopy(self.surface)
        boundary_layer_cell_array = []
        boundary_layer_celltype_array = []
        cell_entity_ids_array = []
        output_surface_cellids_array = []
        quadratic = False
        if quadratic:
            number_layer_points = 2 * grid.n_points
        else:
            number_layer_points = grid.n_points
        output_points = np.zeros((grid.n_points + number_layer_points * self.number_of_sublayers, 3))
        output_points[0:grid.n_points, :] = grid.points.astype(np.float64)
        if self.include_surface_cells:
            for i in range(grid.n_cells):
                #cell = grid.get_cell(i)
                cell = self.surface_cells[i]
                surfacePts = []
                if cell['type'].name == 'TRIANGLE':
                    boundary_layer_celltype_array.append(pv.CellType.TRIANGLE)
                    surfacePts.append(cell['point_ids'][0])
                    surfacePts.append(cell['point_ids'][1])
                    surfacePts.append(cell['point_ids'][2])
                elif cell['type'].name == 'QUAD':
                    boundary_layer_celltype_array.append(pv.CellType.QUAD)
                    surfacePts.append(cell['point_ids'][0])
                    surfacePts.append(cell['point_ids'][1])
                    surfacePts.append(cell['point_ids'][2])
                    surfacePts.append(cell['point_ids'][3])
                elif cell['type'].name == 'QUADRATIC_TRIANGLE':
                    pass
                else:
                    raise ValueError("Unsupported cell type.")
                surfacePts.insert(0, len(surfacePts))
                boundary_layer_cell_array.append(surfacePts)
                cell_entity_ids_array.append(self.inner_surface_cell_entity_id)
                output_surface_cellids_array.append(i)
        relaxation = self.relaxation
        initial_number_of_steps = int(self.number_of_substeps / 100)
        intermediate_number_of_steps = int(self.number_of_substeps / 10)
        final_number_of_steps = int(self.number_of_substeps - initial_number_of_steps - intermediate_number_of_steps)
        self.build_warp_vectors(grid)
        self.increment_warp_vectors(grid, initial_number_of_steps, relaxation)
        tangled_cell_ids = self.check_tangle(grid)
        maximum_iterations = int(self.number_of_substeps / 10)
        iteration = 0
        while len(tangled_cell_ids) > 0 and iteration < maximum_iterations:
            self.increment_warp_vectors(grid, intermediate_number_of_steps, relaxation)
            tangled_cell_ids = self.check_tangle(grid)
            self.local_untangle_fast(grid, tangled_cell_ids, self.local_correction_factor)
            tangled_cell_ids = self.check_tangle(grid)
            iteration += 1
            if len(tangled_cell_ids) == 0:
                print("try last : ", iteration)
                self.increment_warp_vectors(grid, final_number_of_steps, relaxation)
                tangled_cell_ids = self.check_tangle(grid)
        input_points = np.array(grid.points, dtype=np.float64)
        for i in tqdm.trange(self.number_of_sublayers, desc="Generating Boundary Layers", leave=False):
            warp_points = self.init_warp_points(input_points, i, False)
            if np.any(np.linalg.norm(warp_points - input_points, axis=1) == 0.0):
                print("Warp points not initialized correctly.")
            output_points[grid.n_points + i * number_layer_points:grid.n_points + (i + 1) * number_layer_points, :] = warp_points
            for j in range(grid.n_cells):
                #cell = grid.get_cell(j)
                cell = self.surface_cells[j]
                if cell['type'].name == 'TRIANGLE' or cell['type'].name == 'QUAD':
                    prismNPts = len(cell['point_ids']) * 2
                    prismPts = []
                    quadNPts = 4
                    quadPts = []
                    for k in range(len(cell['point_ids'])):
                        prismPts.append(cell['point_ids'][k] + i*number_layer_points)
                    for k in range(len(cell['point_ids'])):
                        prismPts.append(cell['point_ids'][k] + (i + 1) * number_layer_points)
                    prismPts.insert(0, prismNPts)
                    boundary_layer_cell_array.append(prismPts)
                    cell_entity_ids_array.append(self.volume_cell_entity_id)
                    output_surface_cellids_array.append(-1)
                    if cell['type'].name == 'TRIANGLE':
                        boundary_layer_celltype_array.append(pv.CellType.WEDGE)
                    elif cell['type'].name == 'QUAD':
                        boundary_layer_celltype_array.append(pv.CellType.HEXAHEDRON)
                    if self.include_sidewall_cells:
                        pass
                elif cell['type'].name == 'QUADRATIC_TRIANGLE':
                    if self.include_sidewall_cells:
                        pass
                    pass
                else:
                    raise ValueError("Unsupported cell type.")
            if self.include_surface_cells:
                if i == self.number_of_sublayers - 1:
                    for j in range(grid.n_cells):
                        #cell = grid.get_cell(j)
                        cell = self.surface_cells[j]
                        surfacePts = []
                        if cell['type'].name == 'TRIANGLE':
                            boundary_layer_celltype_array.append(pv.CellType.TRIANGLE)
                            surfacePts.append(cell['point_ids'][0] + i * number_layer_points + number_layer_points)
                            surfacePts.append(cell['point_ids'][1] + i * number_layer_points + number_layer_points)
                            surfacePts.append(cell['point_ids'][2] + i * number_layer_points + number_layer_points)
                        elif cell['type'].name == 'QUAD':
                            boundary_layer_celltype_array.append(pv.CellType.QUAD)
                            surfacePts.append(cell['point_ids'][0] + i * number_layer_points + number_layer_points)
                            surfacePts.append(cell['point_ids'][1] + i * number_layer_points + number_layer_points)
                            surfacePts.append(cell['point_ids'][2] + i * number_layer_points + number_layer_points)
                            surfacePts.append(cell['point_ids'][3] + i * number_layer_points + number_layer_points)
                        elif cell['type'].name == 'QUADRATIC_TRIANGLE':
                            pass
                        else:
                            raise ValueError("Unsupported cell type.")
                        surfacePts.insert(0, len(surfacePts))
                        boundary_layer_cell_array.append(surfacePts)
                        cell_entity_ids_array.append(self.outer_surface_cell_entity_id)
                        output_surface_cellids_array.append(j)
        all_cell_array = []
        for cell in boundary_layer_cell_array:
            all_cell_array.extend(cell)
        boundary_layer_cell_array = np.array(all_cell_array).ravel()
        layers = pv.UnstructuredGrid(boundary_layer_cell_array, boundary_layer_celltype_array, output_points)
        layers.cell_data['EntityID'] = np.array(cell_entity_ids_array)
        layers = self.tetra_filter(layers)
        self.layers = layers
        if self.combine:
            combined, interior = self.build_combined(layers, tolerance=self.tolerance)
        else:
            combined = None
            interior = None
        return combined, interior, layers

    def build_warp_vectors(self, grid):
        self.warp_vectors = np.array(grid.point_data[self.warp_vector_array], dtype=np.float64)
        if self.negate_warp_vectors:
            self.warp_vectors = -self.warp_vectors
        layer_thickness = 0.0
        if self.constant_thickness:
            layer_thickness = self.layer_thickness
        elif self.use_warp_vector_magnitude_as_thickness:
            layer_thickness = np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
        else:
            layer_thickness = self.layer_thickness_array * self.layer_thickness_ratio
        if self.max_layer_thickness is not None:
            layer_thickness = np.minimum(layer_thickness, self.max_layer_thickness)
        magnitudes = np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
        if np.any(magnitudes == 0.0):
            raise ValueError("Warp vectors have zero magnitude.")
        self.warp_vectors = (self.warp_vectors / magnitudes) * layer_thickness
        return

    def increment_warp_vectors(self, grid, number_of_substeps, relaxation):
        input_points = np.array(grid.points, dtype=np.float64)
        base_points = np.array(grid.points, dtype=np.float64)
        warp_points = np.zeros_like(base_points)
        for step in range(number_of_substeps):
            warp_points = self.inrcemental_warp_points(grid, base_points, warp_points, step, number_of_substeps, relaxation)
            base_points = deepcopy(warp_points)
        self.warp_vectors = warp_points - input_points
        layer_thickness = np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
        self.warp_vectors = self.warp_vectors / layer_thickness
        self.warp_vectors = self.warp_vectors * layer_thickness
        return

    def inrcemental_warp_points(self, grid, base_points, warp_points, step, number_of_steps, relaxation):
        layer_thickness = np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
        warp_vectors = self.warp_vectors / layer_thickness.reshape(-1, 1)
        layer_thickness /= number_of_steps
        warp_points = base_points + warp_vectors * layer_thickness
        edges = grid.extract_feature_edges(boundary_edges=True, manifold_edges=False,
                                                 feature_edges=False, non_manifold_edges=False)
        edge_pt_ids = []
        all_points = np.array(grid.points, dtype=np.float64)
        for i in range(edges.n_points):
            edge_pt_ids.append(np.argwhere(np.all(edges.points[i, :] == grid.points, axis=1)).flatten()[0])
        for i in range(grid.points.shape[0]):
            #neighbors = grid.point_neighbors(i)
            neighbors = self.point_neighbors[i]
            # Determine if point is on the edge of the surface
            if i in edge_pt_ids:
                continue
            barycenter = all_points[neighbors].mean(axis=0)
            warp_points[i, :] += relaxation * (barycenter - warp_points[i, :])
        return warp_points

    def check_tangle(self, grid):
        found = np.zeros(grid.n_cells)
        all_points = np.array(grid.points, dtype=np.float64)
        #for i in range(grid.n_cells):
        #    #cell = grid.get_cell(i)
        #    cell = self.surface_cells[i]
        #    cell_points = all_points[cell['point_ids'], :]
        #    warped_points = cell_points + self.warp_vectors[cell['point_ids'], :]
        #    base_normal = self.compute_normal(cell_points)
        #    warped_normal = self.compute_normal(warped_points)
        #    base_area = self.compute_area(cell_points)
        #    warped_area = self.compute_area(warped_points)
        #    if np.dot(base_normal, warped_normal) < 0.0 or warped_area / base_area < 0.1:
        #        found.append(i)
        cell_points = all_points[self.surface_cell_array, :]
        warped_points = cell_points + self.warp_vectors[self.surface_cell_array, :]
        base_normals = self.compute_normals(cell_points)
        warped_normals = self.compute_normals(warped_points)
        base_area = self.compute_areas(cell_points)
        warped_area = self.compute_areas(warped_points)
        normal_dots = (np.linalg.norm(base_normals * warped_normals, axis=1) < 0.0).flatten()
        found[normal_dots] = 1
        area_ratios = ((warped_area / base_area) <= 0.1).flatten()
        found[area_ratios] = 1
        found = np.argwhere(found == 1).flatten().tolist()
        return found

    def local_untangle(self, grid, tangled_cell_ids, alpha):
        correction_array = np.zeros((grid.n_points, 3), dtype=np.float64)
        for i in range(grid.n_cells):
            if i not in tangled_cell_ids:
                continue
            #cell = grid.get_cell(i)
            cell = self.surface_cells[i]
            warp_vectors = self.warp_vectors[cell['point_ids'], :]
            # SimVascular Implementation only does 1 iteration?
            for j in range(1):
                w1s = np.zeros(3)
                w2s = np.zeros(3)
                w3s = np.zeros(3)
                neighbors = grid.point_cell_ids(cell['point_ids'][0])
                for k in neighbors:
                    #neighbor_cell = grid.get_cell(k)
                    neighbor_cell = self.surface_cells[k]
                    n = self.compute_normal(neighbor_cell['points'])
                    w1s += warp_vectors[0, :] - np.dot(warp_vectors[0, :], n) * n
                w1s_mag = np.linalg.norm(w1s)
                w1s = w1s / w1s_mag
                neighbors = grid.point_cell_ids(cell['point_ids'][1])
                for k in neighbors:
                    #neighbor_cell = grid.get_cell(k)
                    neighbor_cell = self.surface_cells[k]
                    n = self.compute_normal(neighbor_cell['points'])
                    w2s += warp_vectors[1, :] - np.dot(warp_vectors[1, :], n) * n
                w2s_mag = np.linalg.norm(w2s)
                w2s = w2s / w2s_mag
                neighbors = grid.point_cell_ids(cell['point_ids'][2])
                for k in neighbors:
                    #neighbor_cell = grid.get_cell(k)
                    neighbor_cell = self.surface_cells[k]
                    n = self.compute_normal(neighbor_cell['points'])
                    w3s += warp_vectors[2, :] - np.dot(warp_vectors[2, :], n) * n
                w3s_mag = np.linalg.norm(w3s)
                w3s = w3s / w3s_mag
                correction_array[cell['point_ids'][0], :] += alpha * w2s + alpha * w3s
                correction_array[cell['point_ids'][1], :] += alpha * w1s + alpha * w3s
                correction_array[cell['point_ids'][2], :] += alpha * w1s + alpha * w2s
        layer_thickness = np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
        self.warp_vectors = (self.warp_vectors + correction_array)
        self.warp_vectors = self.warp_vectors / np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
        self.warp_vectors = self.warp_vectors * layer_thickness
        if np.any(layer_thickness == 0.0):
            warnings.warn("Warp vectors have zero magnitude.")
        return

    def local_untangle_fast(self, grid, tangled_cell_ids, alpha):
        correction_array = np.zeros((grid.n_points, 3), dtype=np.float64)
        all_cells = self.surface.faces.reshape(-1, 4)[:, 1:]

        for i in tqdm.trange(grid.n_cells, desc="Local Untangle", leave=False):
            if i not in tangled_cell_ids:
                continue
            #cell = grid.get_cell(i)
            cell = self.surface_cells[i]
            warp_vectors = self.warp_vectors[cell['point_ids'], :]
            # SimVascular Implementation only does 1 iteration?
            for j in range(1):
                w1s = np.zeros(3)
                w2s = np.zeros(3)
                w3s = np.zeros(3)
                #neighbors = grid.point_cell_ids(cell['point_ids'][0])
                neighbors = self.point_cells[cell['point_ids'][0]]
                for k in neighbors:
                    #neighbor_cell = grid.get_cell(k)
                    neighbor_cell = self.surface_cells[k]
                    n = self.compute_normal(neighbor_cell['points'])
                    w1s += warp_vectors[0, :] - np.dot(warp_vectors[0, :], n) * n
                w1s_mag = np.linalg.norm(w1s)
                w1s = w1s / w1s_mag
                #neighbors = grid.point_cell_ids(cell['point_ids'][1])
                neighbors = self.point_cells[cell['point_ids'][1]]
                for k in neighbors:
                    #neighbor_cell = grid.get_cell(k)
                    neighbor_cell = self.surface_cells[k]
                    n = self.compute_normal(neighbor_cell['points'])
                    w2s += warp_vectors[1, :] - np.dot(warp_vectors[1, :], n) * n
                w2s_mag = np.linalg.norm(w2s)
                w2s = w2s / w2s_mag
                #neighbors = grid.point_cell_ids(cell['point_ids'][2])
                neighbors = self.point_cells[cell['point_ids'][2]]
                for k in neighbors:
                    #neighbor_cell = grid.get_cell(k)
                    neighbor_cell = self.surface_cells[k]
                    n = self.compute_normal(neighbor_cell['points'])
                    w3s += warp_vectors[2, :] - np.dot(warp_vectors[2, :], n) * n
                w3s_mag = np.linalg.norm(w3s)
                w3s = w3s / w3s_mag
                correction_array[cell['point_ids'][0], :] += alpha * w2s + alpha * w3s
                correction_array[cell['point_ids'][1], :] += alpha * w1s + alpha * w3s
                correction_array[cell['point_ids'][2], :] += alpha * w1s + alpha * w2s
        layer_thickness = np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
        self.warp_vectors = (self.warp_vectors + correction_array)
        self.warp_vectors = self.warp_vectors / np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
        self.warp_vectors = self.warp_vectors * layer_thickness
        if np.any(layer_thickness == 0.0):
            warnings.warn("Warp vectors have zero magnitude.")
        return

    def init_warp_points(self, input_points, sub_layer_id, quadratic):
        total_layer_zero_sublayer_ratio = 0.0
        for i in range(self.number_of_sublayers):
            total_layer_zero_sublayer_ratio += self.sublayer_ratio ** (self.number_of_sublayers - i - 1)
        sublayer_offset_ratio = 0.0
        for i in range(0, sub_layer_id):
            sublayer_offset_ratio += (self.sublayer_ratio ** (self.number_of_sublayers - i - 1))
        sublayer_offset_ratio /= total_layer_zero_sublayer_ratio
        sublayer_thickness_ratio = (self.sublayer_ratio ** (self.number_of_sublayers - sub_layer_id - 1)) / total_layer_zero_sublayer_ratio
        if not quadratic:
            warped_points = np.zeros((input_points.shape[0], 3))
            layer_thickness = np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
            warp_vectors = self.warp_vectors / layer_thickness
            sublayer_offset = sublayer_offset_ratio * layer_thickness
            sublayer_thickness = sublayer_thickness_ratio * layer_thickness
            warped_points[:, :] = input_points + warp_vectors * (sublayer_offset + sublayer_thickness)
        else:
            warped_points = np.zeros((2 * input_points.shape[0], 3))
            layer_thickness = np.linalg.norm(self.warp_vectors, axis=1).reshape(-1, 1)
            warp_vectors = self.warp_vectors / layer_thickness
            sublayer_offset = sublayer_offset_ratio * layer_thickness
            sublayer_thickness = sublayer_thickness_ratio * layer_thickness
            warped_points[0:input_points.shape[0], :] = input_points + 0.5 * warp_vectors * (sublayer_offset + sublayer_thickness)
            warped_points[input_points.shape[0]:, :] = input_points + warp_vectors * (sublayer_offset + sublayer_thickness)
        return warped_points

    def compute_normal(self, pts):
        n = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        return n / np.linalg.norm(n)

    def compute_area(self, pts):
        n = np.cross(pts[1] - pts[0], pts[2] - pts[0])
        return 0.5 * np.linalg.norm(n)

    def compute_normals(self, pts):
        n = np.zeros((pts.shape[0], 3), dtype=np.float64)
        n[:, :] = np.cross(pts[:, 1, :] - pts[:, 0, :], pts[:, 2, :] - pts[:, 0, :])
        return n

    def compute_areas(self, pts):
        n = np.zeros((pts.shape[0], 3), dtype=np.float64)
        n[:, :] = np.cross(pts[:, 1, :] - pts[:, 0, :], pts[:, 2, :] - pts[:, 0, :])
        return 0.5 * np.linalg.norm(n, axis=1)

    def build_connections(self, layers):
        connections = {'vol_to_vol':{}, 'vol_to_faces':{}, 'faces_to_vol':{}, "faces_to_faces":{}}
        volume_cells = np.argwhere(layers.cell_data['EntityID'] == self.volume_cell_entity_id).flatten()
        volume_cells = layers.extract_cells(volume_cells)

    def tetra_filter(self, layers):
        """
        This function filters the boundary layer mesh to yield a tetrahedral mesh.
        We assume that the volume elements of the unstructured grid are well-posed
        in that they define a closed volume and represent a valid element template.
        However, we do not make assumptions regarding orientation of the volume
        elements.
        :param layers:
        :return:
        """
        tetra_ordering = list(permutations([1, 2, 3, 4]))
        inner_submesh_cells = []
        inner_submesh_cells_type = []
        inner_submesh_cells_entity = []
        outer_submesh_cells = []
        outer_submesh_cells_type = []
        outer_submesh_cells_entity = []
        inner_submesh_ids = np.argwhere(layers.cell_data['EntityID'] == self.inner_surface_cell_entity_id).flatten()
        outer_submesh_ids = np.argwhere(layers.cell_data['EntityID'] == self.outer_surface_cell_entity_id).flatten()
        for i in inner_submesh_ids:
            cell = layers.get_cell(i)
            nodes = cell.point_ids
            nodes.insert(0, len(nodes))
            inner_submesh_cells.append(nodes)
            inner_submesh_cells_type.append(cell.type)
            inner_submesh_cells_entity.append(self.inner_surface_cell_entity_id)
        for i in outer_submesh_ids:
            cell = layers.get_cell(i)
            nodes = cell.point_ids
            nodes.insert(0, len(nodes))
            outer_submesh_cells.append(nodes)
            outer_submesh_cells_type.append(cell.type)
            outer_submesh_cells_entity.append(self.outer_surface_cell_entity_id)
        volume_submesh = layers.extract_cells(np.argwhere(layers.cell_data['EntityID'] == self.volume_cell_entity_id).flatten())
        volume_submesh = volume_submesh.triangulate(progress_bar=True)
        volume_submesh_quality = compute_cell_quality(volume_submesh).cell_data['CellQuality'].flatten()
        #if np.any(volume_submesh_quality < 0):
        #    old_tetra_cells = deepcopy(volume_submesh.cells.reshape(-1, 5))
        #    new_tetra_cells = deepcopy(volume_submesh.cells.reshape(-1, 5))
        #    tetra_cell_type = [pv.CellType.TETRA for i in range(new_tetra_cells.shape[0])]
        #    output_points = np.array(volume_submesh.points)
        #    for i in tqdm.trange(len(tetra_ordering), desc="Reordering Tetrahedra", leave=False):
        #        ordering = tetra_ordering[i]
        #        ids = np.argwhere(volume_submesh_quality < 0)
        #        #changed = tetra_cells[ids, ordering].copy()
        #        new_tetra_cells[ids.flatten(), 1:] = old_tetra_cells[ids, ordering]
        #        new_volume_submesh = pv.UnstructuredGrid(new_tetra_cells.ravel(), tetra_cell_type, output_points)
        #        quality = new_volume_submesh.compute_cell_quality().cell_data['CellQuality'].flatten()
        #        if not np.any(quality < 0):
        #            break
        volume_submesh_cells = volume_submesh.cells.reshape(-1, 5).tolist()
        volume_submesh_cells_type = [pv.CellType.TETRA] * len(volume_submesh_cells)
        volume_submesh_entity = [self.volume_cell_entity_id] * len(volume_submesh_cells)
        all_cells = []
        all_cells.extend(inner_submesh_cells)
        all_cells.extend(volume_submesh_cells)
        all_cells.extend(outer_submesh_cells)
        all_cells_type = []
        all_cells_type.extend(inner_submesh_cells_type)
        all_cells_type.extend(volume_submesh_cells_type)
        all_cells_type.extend(outer_submesh_cells_type)
        all_cells_entity = []
        all_cells_entity.extend(inner_submesh_cells_entity)
        all_cells_entity.extend(volume_submesh_entity)
        all_cells_entity.extend(outer_submesh_cells_entity)
        tmp = []
        for i in range(len(all_cells)):
            tmp.extend(all_cells[i])
        all_cells = np.array(tmp).ravel()
        all_cells_type = np.array(all_cells_type)
        all_cells_entity = np.array(all_cells_entity)
        all_points = np.array(layers.points)
        tetra_mesh = pv.UnstructuredGrid(all_cells, all_cells_type, all_points)
        tetra_mesh.cell_data['EntityID'] = all_cells_entity
        return tetra_mesh

    def build_combined(self, layers, tolerance=1E-14):
        # The points of the mesh surface need to be shifted to correct the mesh for the
        # inclusion of the boundary layers
        layers = deepcopy(layers)
        inner_surface_cells = np.argwhere(layers.cell_data['EntityID'] == self.inner_surface_cell_entity_id).flatten()
        outer_surface_cells = np.argwhere(layers.cell_data['EntityID'] == self.outer_surface_cell_entity_id).flatten()
        inner_surface_pt_ids = []
        outer_surface_pt_ids = []
        for i in tqdm.trange(len(inner_surface_cells), desc="Collecting Inner Surface Points", leave=False):
            cell_id = inner_surface_cells[i]
            cell = layers.get_cell(cell_id)
            inner_surface_pt_ids.extend(cell.point_ids)
        for i in tqdm.trange(len(outer_surface_cells), desc="Collecting Outer Surface Points", leave=False):
            cell_id = outer_surface_cells[i]
            cell = layers.get_cell(cell_id)
            outer_surface_pt_ids.extend(cell.point_ids)
        outer = layers.extract_cells(outer_surface_cells)
        layers_tree = cKDTree(layers.points)
        dists, inds = layers_tree.query(outer.points, k=2)
        if np.any(dists[:, 1] == 0.0):
            warnings.warn("Internal Layer Points coincident with Outer Surface Points.")
        original_point_ids = outer.point_data['vtkOriginalPointIds']
        outer = outer.extract_surface()
        outer_tree = cKDTree(outer.points)

        boundaries = outer.extract_feature_edges(boundary_edges=True, manifold_edges=False, feature_edges=False,
                                                 non_manifold_edges=False)
        boundaries = boundaries.split_bodies()
        caps = []
        for i in range(len(boundaries)):
            cap = remesh_surface_2d(boundaries[i])
            #_, outer_cap_inds = outer_tree.query(cap.points)
            #cap.points[:boundaries[i].n_points, :] = outer.points[outer_cap_inds[:boundaries[i].n_points], :]
            caps.append(cap)
        caps.insert(0, outer)
        outer_total = pv.merge(caps, merge_points=True, progress_bar=True)
        outer_total = outer_total.clean(tolerance=tolerance, progress_bar=True)
        if not outer_total.is_manifold:
            raise ValueError("Merged Outer Surface is not manifold.")
        else:
            print("Merged Outer Surface is manifold.")
        outer_total_quality = compute_cell_quality(outer_total).cell_data["CellQuality"]
        if np.any(outer_total_quality < 0.0):
            warnings.warn("Outer Surface has inverted triangles.")
        interior = tetgen.TetGen(outer_total)
        interior.tetrahedralize(switches='pq1.1/10MVYSJ')
        mesh_interior = interior.grid
        if self.remesh_vol:
            print("Remeshing interior volume...")
            mesh_interior = remesh_volume(mesh_interior, nosurf=True)
            print("Remeshing interior complete.")
        mesh_interior = mesh_interior.triangulate().clean()
        mesh_interior_surface = mesh_interior.extract_surface()
        if not mesh_interior_surface.is_manifold:
            warnings.warn("Interior Surface is not manifold.")
        else:
            print("Interior Surface is manifold.")
        layers_volume = layers.extract_cells(np.argwhere(layers.cell_data['EntityID'] == self.volume_cell_entity_id).flatten())
        layers_volume = layers_volume.triangulate().clean()
        layers_surface = layers_volume.extract_surface()
        if not layers_surface.is_manifold:
            warnings.warn("Layers Surface is not manifold.")
        else:
            print("Layers Surface is manifold.")
        # Check Quality of Combined Mesh
        quality = compute_cell_quality(mesh_interior).cell_data['CellQuality'].flatten()
        if np.any(quality < 0):
            warnings.warn("Combined Mesh has negative jacobian at {}".format(np.argwhere(quality < 0).flatten()))
            print("Attempting to fix negative jacobians via volume meshing.")
        """
        layer_to_outer_dists, layer_to_outer_inds = outer_tree.query(layers_volume.points)
        outer_to_layer = {}
        tolerance = 1E-6
        for i in tqdm.trange(len(layer_to_outer_inds), desc="Mapping Layer Points to Outer Surface", leave=False):
            if layer_to_outer_dists[i] < tolerance:
                if outer_to_layer.get(layer_to_outer_inds[i], None) is None:
                    outer_to_layer[layer_to_outer_inds[i]] = i
                else:
                    warnings.warn("Multiple Layer Points Mapped to Single Outer Point.")
        interior_to_outer_dists, interior_to_outer_inds = outer_tree.query(mesh_interior.points)
        outer_to_interior = {}
        for i in tqdm.trange(len(interior_to_outer_inds), desc="Mapping Interior Points to Outer Surface", leave=False):
            if interior_to_outer_dists[i] < tolerance:
                if outer_to_interior.get(interior_to_outer_inds[i], None) is None:
                    outer_to_interior[interior_to_outer_inds[i]] = i
                else:
                    warnings.warn("Multiple Interior Points Mapped to Single Outer Point.")
        layer_cells = layers_volume.cells.reshape(-1, 5)[:, 1:]
        """

        combined_mesh = pv.merge([layers_volume, mesh_interior], merge_points=True, tolerance=tolerance)
        #combined_mesh = combined_mesh.clean()
        #combined_mesh = pv.merge([layers_volume, mesh_interior])
        #combined_mesh = combined_mesh.triangulate().clean()
        combined_mesh_surface = combined_mesh.extract_surface()
        if not combined_mesh_surface.is_manifold:
            non_manifold_edges = combined_mesh_surface.extract_feature_edges(boundary_edges=False, manifold_edges=False,
                                                                             feature_edges=False, non_manifold_edges=True)
            warnings.warn("Combined Mesh Surface is not manifold.\n {} Non-Manifold Edges at: {}".format(non_manifold_edges.n_cells, non_manifold_edges.cell_data["vtkOriginalCellIds"]) )
        else:
            print("Combined Mesh Surface is manifold.")
        # combined_mesh = pv.merge([layers_volume, mesh_interior], merge_points=True)
        # Replace the pv.merge function with a custom merge function
        return combined_mesh, mesh_interior

def append_boundary_layers(mesh, layers, inner_surface_cell_entity_id=1, outer_surface_cell_entity_id=2,
                           compute_inner_outer_mapping=False, scale=20):
    # The points of the mesh surface need to be shifted to correct the mesh for the
    # inclusion of the boundary layers
    mesh = deepcopy(mesh)
    layers = deepcopy(layers)
    inner_surface_cells = np.argwhere(layers.cell_data['EntityID'] == inner_surface_cell_entity_id).flatten()
    outer_surface_cells = np.argwhere(layers.cell_data['EntityID'] == outer_surface_cell_entity_id).flatten()
    inner_surface_pt_ids = []
    outer_surface_pt_ids = []
    for i in tqdm.trange(len(inner_surface_cells), desc="Collecting Inner Surface Points", leave=False):
        cell_id = inner_surface_cells[i]
        cell = layers.get_cell(cell_id)
        inner_surface_pt_ids.extend(cell.point_ids)
    for i in tqdm.trange(len(outer_surface_cells), desc="Collecting Outer Surface Points", leave=False):
        cell_id = outer_surface_cells[i]
        cell = layers.get_cell(cell_id)
        outer_surface_pt_ids.extend(cell.point_ids)
    inner_surface_pt_ids = np.unique(np.array(inner_surface_pt_ids).flatten())
    outer_surface_pt_ids = np.unique(np.array(outer_surface_pt_ids).flatten())
    inner_tree = cKDTree(layers.points[inner_surface_pt_ids])
    outer_tree = cKDTree(layers.points[outer_surface_pt_ids])
    mesh_tree = cKDTree(mesh.points)
    if compute_inner_outer_mapping:
        # Map the inner surface points to the outer surface points
        cell_volumes = mesh.compute_cell_sizes().cell_data['Volume']
        length = (cell_volumes.mean()**(1/3)) * scale
        sparse = inner_tree.sparse_distance_matrix(outer_tree, length, output_type='coo_matrix')
        row, col = min_weight_full_bipartite_matching(sparse)
        inner_to_outer_idx = outer_surface_pt_ids[col]
    else:
        inner_to_outer_idx = deepcopy(outer_surface_pt_ids)
    # Check that inner <-> outer mapping is 1:1
    _, counts = np.unique(inner_to_outer_idx, return_counts=True)
    if np.any(counts > 1):
        raise ValueError("Inner to Outer Surface Point Mapping is not 1:1.")
    # Map the inner surface points to the mesh surface points
    inner_to_mesh_dist, inner_to_mesh_idx = mesh_tree.query(layers.points[inner_surface_pt_ids])
    # Check that inner <-> mesh mapping is 1:1
    _, _, counts = np.unique(inner_to_mesh_idx, return_inverse=True, return_counts=True)
    if np.any(counts > 1):
        raise ValueError("Inner to Mesh Surface Point Mapping is not 1:1.")
    # Check that inner <-> mesh mapping is conforming
    if np.any(~np.isclose(inner_to_mesh_dist, 0.0)):
        raise ValueError("Inner to Mesh Surface Point Mapping is not conforming.\n Could result in mesh corruption.")
    outer = layers.extract_cells(outer_surface_cells).extract_surface()
    outer_tree = cKDTree(outer.points)
    boundaries = outer.extract_feature_edges(boundary_edges=True,manifold_edges=False,feature_edges=False,non_manifold_edges=False)
    boundaries = boundaries.split_bodies()
    caps = []
    for i in range(len(boundaries)):
        cap = remesh_surface_2d(boundaries[i])
        _, outer_cap_inds = outer_tree.query(cap.points)
        cap.points[:boundaries[i].n_points, :] = outer.points[outer_cap_inds[:boundaries[i].n_points], :]
        caps.append(cap)
    caps.insert(0, outer)
    outer_total = pv.merge(caps)
    interior = tetgen.TetGen(outer_total)
    interior.tetrahedralize()
    mesh = interior.grid
    #mesh.points[inner_to_mesh_idx[inner_surface_pt_ids]] = layers.points[inner_to_outer_idx]
    # Check Quality of Combined Mesh
    quality = compute_cell_quality(mesh).cell_data['CellQuality'].flatten()
    if np.any(quality < 0):
        warnings.warn("Combined Mesh has negative jacobian at {}".format(np.argwhere(quality < 0).flatten()))
        print("Attempting to fix negative jacobians via volume meshing.")
    #combined_mesh = pv.merge([layers, mesh], merge_points=True)
    combined_mesh = pv.merge([layers, mesh], merge_points=True)
    return combined_mesh, mesh
