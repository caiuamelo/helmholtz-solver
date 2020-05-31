import numpy as np
from function_builders import RegionFunctionBuilder
import globals as g
import _fwi_ls


def build_mesh(size_x, size_y, nx, ny):
    build_connectivity_list = g.choose_impl(_fwi_ls.build_connectivity_list, _build_connectivity_list)

    x = np.linspace(0, size_x, nx + 1)
    y = np.linspace(0, size_y, ny + 1)
    xv, yv = np.meshgrid(x, y)
    points = np.vstack([xv.ravel(), yv.ravel()]).T
    connectivity_list = build_connectivity_list(nx, ny)
    mesh = Mesh(nx, ny, points, connectivity_list, size_x, size_y)

    mu = lambda c: 1.0 / np.power(c, 2)
    mu_f = (
        RegionFunctionBuilder(overall_value=mu(c=300.0))
        .set_rectangle_interval_value(
            start_point=(60.0, 60.0), end_point=(140.0, 140.0), value=mu(c=250.0),
        )
        .build()
    )

    eta_f = (
        RegionFunctionBuilder(overall_value=4e-3)
        .set_rectangle_interval_value(
            start_point=(20.0, 20.0), end_point=(180.0, 180.0), value=0.0
        )
        .build()
    )

    mesh.add_source("source 1", (30.0, 30.0), 1.0)
    mesh.set_mu(mu_f)
    mesh.set_eta(eta_f)

    return mesh


def _build_connectivity_list(nx, ny):
    from_ij_to_flat_index = lambda i, j: (nx + 1) * i + j
    connectivity_list = [None] * (nx * ny)
    for i in range(ny):
        for j in range(nx):
            connectivity_list[i * nx + j] = [
                from_ij_to_flat_index(i, j),
                from_ij_to_flat_index(i, j + 1),
                from_ij_to_flat_index(i + 1, j + 1),
                from_ij_to_flat_index(i + 1, j),
            ]
    return np.array(connectivity_list)


class Mesh:
    class Source:
        def __init__(self, name, position, closest_point_id, value):
            self.name = name
            self.position = position
            self.closest_node_id = closest_point_id
            self.value = value

    def __init__(self, nx, ny, points, connectivity_list, size_x, size_y):
        self.points = points
        self.connectivity_list = connectivity_list
        self.n_points = (nx + 1) * (ny + 1)
        self.nx = nx
        self.ny = ny
        self.size_x = size_x
        self.size_y = size_y
        self.points_in_elements = self._build_points_in_elements()
        self._active_sources = []
        self._sources = {}
        self._source_at_element = {}

    def closest_node_id(self, desired_position):
        closest_distance = +np.inf
        closest_id = None
        for i, point in enumerate(self.points):
            d = np.linalg.norm(point - desired_position, ord=2)
            if d < closest_distance:
                closest_distance = d
                closest_id = i
        return closest_id

    def set_mu(self, mu_f):
        self.mu = np.empty(shape=(len(self.points_in_elements),))
        for i, points in enumerate(self.points_in_elements):
            x, y = np.mean(points, axis=0)
            self.mu[i] = mu_f(x, y)

    def set_eta(self, mu_f):
        self.eta = np.empty(shape=(len(self.points_in_elements),))
        for i, points in enumerate(self.points_in_elements):
            x, y = np.mean(points, axis=0)
            self.eta[i] = mu_f(x, y)

    def add_source(self, source_name, position, value):
        closest_point_id = self.closest_node_id(position)
        eid = self._element_id_in_position(position)
        s = Mesh.Source(source_name, position, closest_point_id, value,)
        self._sources[source_name] = s
        self._source_at_element[eid] = s.value
        self._active_sources.append(source_name)

    def source_at_element(self, eid):
        return self._source_at_element.get(eid, 0.0)

    def _element_id_in_position(self, position):
        x, y = position
        for i, element_points in enumerate(self.points_in_elements):
            xini, yini = element_points[0]
            xend, yend = element_points[2]
            if xini <= x <= xend and yini <= y <= yend:
                return i
        assert False, f"No element found in position {position}"

    def _build_points_in_elements(self):
        return np.array(
            [[self.points[pid] for pid in e] for e in self.connectivity_list]
        )
