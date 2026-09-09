"""Duffy-only quadrature for the scalar nearfield notebook comparisons.

This mirrors ProjectPointToReference/GetIntegrationRule in bem/potentialcf.cpp:
five constrained Gauss-Newton steps, a 1e-12 reference-step tolerance, and a
triangle fan around the projection. NGSolve integrates the original curved
surface integrand. No tangent subtraction or analytical correction is added.
"""

from functools import lru_cache

import numpy as np
import ngsolve as ng
from ngsolve.fem import IntegrationRule


_corners = {
    ng.TRIG: np.array(((0., 0.), (1., 0.), (0., 1.))),
    ng.QUAD: np.array(((0., 0.), (1., 0.), (1., 1.), (0., 1.))),
}


def _minimize_on_reference(a, b, element_type):
    """Minimize 0.5 * uv @ a @ uv + b @ uv on the reference element."""
    if a[0, 0] > 0 and np.linalg.det(a) > 0:
        uv = np.linalg.solve(a, -b)
        inside = (np.all(uv > 0) and
                  (uv.sum() < 1 if element_type == ng.TRIG else np.all(uv < 1)))
        if inside:
            return uv

    candidates = []
    corners = _corners[element_type]
    for start, end in zip(corners, np.roll(corners, -1, axis=0)):
        direction = end - start
        curvature = direction @ a @ direction
        slope = (a @ start + b) @ direction
        t = np.clip(-slope / curvature, 0., 1.) if curvature > 0 else float(
            0.5 * curvature + slope <= 0)
        uv = start + t * direction
        candidates.append((0.5 * uv @ a @ uv + b @ uv, uv))
    return min(candidates, key=lambda candidate: candidate[0])[1]


def project_point_to_reference(point, trafo, element_type):
    """Local constrained projection, with the same stopping rule as C++."""
    uv = _corners[element_type].mean(axis=0)
    for _ in range(5):
        mip = trafo(*uv)
        jac = np.array(mip.jacobi)
        a = jac.T @ jac
        b = -jac.T @ (point - np.array(mip.point) + jac @ uv)
        candidate = _minimize_on_reference(a, b, element_type)
        step = np.linalg.norm(candidate - uv)
        uv = candidate
        if step <= 1e-12:
            break
    return uv


def duffy_rule(element_type, projection, order):
    """Build the three-sector triangle or four-sector quadrilateral rule."""
    segment = IntegrationRule(ng.SEGM, order)
    nodes = np.array([ip.point[0] for ip in segment])
    weights = np.array(segment.weights)
    s, t = np.meshgrid(nodes, nodes, indexing='ij')
    tensor_weights = weights[:, None] * weights[None, :] * (1 - t)
    points, mapped_weights = [], []
    corners = _corners[element_type]
    for start, end in zip(corners, np.roll(corners, -1, axis=0)):
        sides = np.column_stack((start - projection, end - projection))
        determinant = np.linalg.det(sides)
        if abs(determinant) <= 1e-12:
            continue
        uv = (start + (s * (1 - t))[..., None] * (end - start)
              + t[..., None] * (projection - start))
        points.extend(map(tuple, uv.reshape(-1, 2)))
        mapped_weights.extend((tensor_weights * determinant).ravel())
    return IntegrationRule(points=points, weights=mapped_weights)


class DuffyQuadrature:
    """Integrate scalar CFs with Duffy near the target and Gauss farther away.

    Create one instance for a fixed source mesh/region and reuse it for the
    distance and order sweeps. Recreate it after changing the mesh geometry.
    The target is a physical coordinate triple; no target mesh is required.
    The nearfield experiments use off-surface targets, not boundary traces.
    """

    def __init__(self, mesh, definedon=None):
        self.mesh = mesh
        self.region = mesh.Boundaries('.*') if definedon is None else definedon
        if mesh.dim != 3 or self.region.VB() != ng.BND:
            raise ValueError('DuffyQuadrature requires a boundary region in 3D')
        self.nelements = len(list(mesh.Elements(ng.BND)))
        self.elements = list(self.region.Elements())
        self.trafos = [mesh.GetTrafo(element) for element in self.elements]
        centers, sizes = [], []
        for element, trafo in zip(self.elements, self.trafos):
            if element.type not in _corners:
                raise ValueError('DuffyQuadrature supports triangles and quads')
            mip = trafo(*_corners[element.type].mean(axis=0))
            centers.append(np.array(mip.point))
            sizes.append(np.linalg.norm(np.array(mip.jacobi)))
        self.centers = np.array(centers).reshape(-1, 3)
        self.sizes = np.array(sizes)

    @lru_cache(maxsize=16)
    def _partition(self, point):
        """Cache only target-dependent projections and element masks."""
        point = np.array(point)
        is_near = np.linalg.norm(self.centers - point, axis=1) < self.sizes
        far = ng.BitArray(self.nelements)
        far.Clear()
        near = []
        for element, trafo, selected in zip(self.elements, self.trafos, is_near):
            if selected:
                mask = ng.BitArray(self.nelements)
                mask.Clear()
                mask[element.nr] = True
                projection = project_point_to_reference(point, trafo, element.type)
                near.append((element.type, projection, mask))
            else:
                far[element.nr] = True
        return far, near

    def integrate(self, integrand, point, order):
        """Integrate an SL/DL integrand using the absolute quadrature order."""
        point = np.asarray(point, dtype=float)
        if point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError('Expected one point with three finite coordinates')
        far, near = self._partition(tuple(point))
        rules = {et: IntegrationRule(et, order) for et in _corners}
        value = ng.Integrate(integrand * ng.ds(
            definedon=self.region, definedonelements=far, intrules=rules), self.mesh)
        for element_type, projection, mask in near:
            rule = duffy_rule(element_type, projection, order)
            value += ng.Integrate(integrand * ng.ds(
                definedon=self.region, definedonelements=mask,
                intrules={element_type: rule}), self.mesh)
        return value
