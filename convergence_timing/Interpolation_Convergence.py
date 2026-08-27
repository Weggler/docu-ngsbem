#!/usr/bin/env python
# coding: utf-8

"""Generate convergence tables for projection-based interpolation on the unit sphere.

The script runs two h-convergence experiments:

* Pi^1 / Pi^0 for scalar traces in H^{1/2}(Gamma) and H^{-1/2}(Gamma)
* Pi^d for tangential traces in H^{-1/2}(div_Gamma, Gamma)

The results are written to

* pi1_pi0_results.csv
* pid_results.csv

The geometry is refined uniformly. For every polynomial order p, the nominal mesh
size is h_l = initial_maxh / 2**l and the curved surface representation is made
isoparametric by ``mesh.Curve(p)`` on every refinement level.
"""

from csv import DictWriter
from math import log, pi, sqrt as math_sqrt
from pathlib import Path

import netgen.meshing as meshing
from netgen.occ import OCCGeometry, Sphere
from ngsolve import (
    BND,
    CF,
    Compress,
    Cross,
    GridFunction,
    H1,
    HDivSurface,
    Integrate,
    Mesh,
    Norm,
    SurfaceL2,
    TaskManager,
    exp,
    specialcf,
    sqrt,
    x,
    y,
    z,
)

PI1_PI0_CSV = Path("pi1_pi0_results.csv")
PID_CSV = Path("pid_results.csv")

ORDERS = range(1, 5)
REFINEMENT_LEVELS = range(5)
INITIAL_MAXH = 0.5
KAPPA = 5.0

SPHERE = Sphere((0, 0, 0), 1)


def _write_results(rows, path):
    """Write a list of dictionaries to CSV, replacing an existing file."""
    if not rows:
        raise ValueError("Cannot write an empty convergence table")

    path = Path(path)
    with path.open("w", newline="") as file:
        writer = DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _new_surface_mesh(initial_maxh):
    """Create the initial surface mesh of the unit sphere."""
    return Mesh(
        OCCGeometry(SPHERE).GenerateMesh(
            maxh=initial_maxh,
            perfstepsend=meshing.MeshingStep.MESHSURFACE,
        )
    )


def _convergence_data(previous, h, errors):
    """Return convergence factors and observed orders for a set of errors."""
    if previous is None:
        nan = float("nan")
        return {
            name: (nan, nan)
            for name in errors
        }

    h_ratio = previous["h"] / h
    return {
        name: (
            previous[name] / value,
            log(previous[name] / value) / log(h_ratio),
        )
        for name, value in errors.items()
    }


def run_pi1_pi0_convergence(
    orders=ORDERS,
    refinement_levels=REFINEMENT_LEVELS,
    initial_maxh=INITIAL_MAXH,
):
    """Run the scalar Pi^1/Pi^0 h-convergence experiment."""
    results = []
    refinement_levels = tuple(refinement_levels)

    for p in orders:
        print(f"Pi1/Pi0: p={p}")
        mesh = _new_surface_mesh(initial_maxh)
        previous = None

        for level in refinement_levels:
            h = initial_maxh / 2**level
            mesh.Curve(p)  # isoparametric geometry approximation
            boundary = mesh.Boundaries(".*")

            # Exact density on Gamma, pulled back to Gamma_h by radial projection.
            radius = sqrt(x * x + y * y + z * z)
            xs = x / radius
            ys = y / radius
            zs = z / radius
            u_exa_pullback = (
                (xs * xs + ys * ys + zs)
                * (xs - xs * xs)
                * (ys - ys * ys)
            )

            # Smooth extension evaluated directly on Gamma_h.
            u_h = (
                (x * x + y * y + z)
                * (x - x * x)
                * (y - y * y)
            )

            fes_h1 = H1(mesh, order=p, definedon=boundary)
            pi1 = GridFunction(fes_h1)

            fes_l2 = SurfaceL2(
                mesh,
                order=p - 1,
                dual_mapping=False,
                definedon=boundary,
            )
            pi0 = GridFunction(fes_l2)

            intorder = 2 * p + 8
            with TaskManager():
                pi1.Interpolate(u_h, definedon=boundary)
                # Pi^0 is the element-local L2 projection.
                pi0.Set(u_h, BND, definedon=boundary, dual=False)

                u_exa_pullback_l2_sq = Integrate(
                    u_exa_pullback**2, mesh, BND, order=intorder
                )
                u_h_l2_sq = Integrate(u_h**2, mesh, BND, order=intorder)

                # Radial geometry error of Gamma_h against the unit sphere.
                err_geo = math_sqrt(
                    Integrate(
                        (radius - 1) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / (4 * pi)
                )

                # Pure FE/PBI error on Gamma_h.
                err_int_1 = math_sqrt(
                    Integrate(
                        (pi1 - u_h) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / u_h_l2_sq
                )
                err_int_0 = math_sqrt(
                    Integrate(
                        (pi0 - u_h) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / u_h_l2_sq
                )

                # Geometry-induced density discrepancy.
                err_geo_density = math_sqrt(
                    Integrate(
                        (u_h - u_exa_pullback) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / u_exa_pullback_l2_sq
                )

                # Combined interpolation + geometry error.
                err_tot_1 = math_sqrt(
                    Integrate(
                        (pi1 - u_exa_pullback) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / u_exa_pullback_l2_sq
                )
                err_tot_0 = math_sqrt(
                    Integrate(
                        (pi0 - u_exa_pullback) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / u_exa_pullback_l2_sq
                )

            errors = {
                "err_geo": err_geo,
                "err_geo_density": err_geo_density,
                "err_int_1": err_int_1,
                "err_tot_1": err_tot_1,
                "err_int_0": err_int_0,
                "err_tot_0": err_tot_0,
            }
            rates = _convergence_data(previous, h, errors)

            row = {
                "p": p,
                "level": level,
                "h": h,
                "N": mesh.GetNE(BND),
                "ndof_Pi1": fes_h1.ndof,
                "ndof_Pi0": fes_l2.ndof,
            }
            for name, value in errors.items():
                cf, alpha = rates[name]
                row[name] = value
                row[f"CF_{name.removeprefix('err_')}"] = cf
                row[f"alpha_{name.removeprefix('err_')}"] = alpha

            results.append(row)
            previous = row

            if level < refinement_levels[-1]:
                mesh.ngmesh.Refine()

    return results


def run_pid_convergence(
    orders=ORDERS,
    refinement_levels=REFINEMENT_LEVELS,
    initial_maxh=INITIAL_MAXH,
    kappa=KAPPA,
):
    """Run the tangential Pi^d h-convergence experiment."""
    results = []
    refinement_levels = tuple(refinement_levels)

    for p in orders:
        print(f"Pi^d: p={p}")
        mesh = _new_surface_mesh(initial_maxh)
        previous = None

        for level in refinement_levels:
            h = initial_maxh / 2**level
            mesh.Curve(p)
            boundary = mesh.Boundaries(".*")

            # Exact density on Gamma, pulled back to Gamma_h.
            radius = sqrt(x * x + y * y + z * z)
            xs = x / radius
            ys = y / radius
            zs = z / radius

            normal_exa_pullback = CF((xs, ys, zs))
            incident_exa_pullback = CF((1, 0, 0)) * exp(-1j * kappa * zs)
            m_exa_pullback = -Cross(
                normal_exa_pullback,
                incident_exa_pullback,
            )

            # Corresponding tangential density intrinsically defined on Gamma_h.
            normal_h = specialcf.normal(3)
            incident_h = CF((1, 0, 0)) * exp(-1j * kappa * z)
            m_h = -Cross(normal_h, incident_h)

            fes_hdiv = Compress(
                HDivSurface(mesh, order=p - 1, complex=True)
            )
            pi_d = GridFunction(fes_hdiv)

            intorder = 2 * p + 8
            with TaskManager():
                pi_d.Set(m_h, BND, definedon=boundary, dual=True)

                m_h_l2_sq = Integrate(
                    Norm(m_h) ** 2,
                    mesh,
                    BND,
                    order=intorder,
                )
                m_exa_pullback_l2_sq = Integrate(
                    Norm(m_exa_pullback) ** 2,
                    mesh,
                    BND,
                    order=intorder,
                )

                # Pure FE/PBI error on Gamma_h.
                err_int_d = math_sqrt(
                    Integrate(
                        Norm(pi_d - m_h) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / m_h_l2_sq
                )

                # Geometry-induced discrepancy of the tangential density.
                err_geo_d = math_sqrt(
                    Integrate(
                        Norm(m_h - m_exa_pullback) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / m_exa_pullback_l2_sq
                )

                # Combined interpolation + geometry error.
                err_tot_d = math_sqrt(
                    Integrate(
                        Norm(pi_d - m_exa_pullback) ** 2,
                        mesh,
                        BND,
                        order=intorder,
                    )
                    / m_exa_pullback_l2_sq
                )

            errors = {
                "err_int_d": err_int_d,
                "err_geo_d": err_geo_d,
                "err_tot_d": err_tot_d,
            }
            rates = _convergence_data(previous, h, errors)

            row = {
                "p": p,
                "level": level,
                "h": h,
                "N": mesh.GetNE(BND),
                "ndof": fes_hdiv.ndof,
            }
            for name, value in errors.items():
                cf, alpha = rates[name]
                row[name] = value
                row[f"CF_{name.removeprefix('err_')}"] = cf
                row[f"alpha_{name.removeprefix('err_')}"] = alpha

            results.append(row)
            previous = row

            if level < refinement_levels[-1]:
                mesh.ngmesh.Refine()

    return results



def run_pi1(order, refinement):
    """Run one Pi^1 convergence sample for pytest.

    ``refinement`` is the number of uniform h-refinement steps applied to the
    same initial surface mesh.  The returned dictionary intentionally follows
    the common convergence-test interface ``{"err": ..., "ndof": ...}``.
    """
    rows = run_pi1_pi0_convergence(
        orders=(order,),
        refinement_levels=range(refinement + 1),
    )
    row = rows[-1]
    return {
        "err": row["err_int_1"],
        "ndof": row["ndof_Pi1"],
    }


def run_pi0(order, refinement):
    """Run one Pi^0 convergence sample for pytest."""
    rows = run_pi1_pi0_convergence(
        orders=(order,),
        refinement_levels=range(refinement + 1),
    )
    row = rows[-1]
    return {
        "err": row["err_int_0"],
        "ndof": row["ndof_Pi0"],
    }


def run_pid(order, refinement):
    """Run one Pi^d convergence sample for pytest."""
    rows = run_pid_convergence(
        orders=(order,),
        refinement_levels=range(refinement + 1),
    )
    row = rows[-1]
    return {
        "err": row["err_int_d"],
        "ndof": row["ndof"],
    }

def main():
    pi1_pi0_results = run_pi1_pi0_convergence()
    _write_results(pi1_pi0_results, PI1_PI0_CSV)
    print(f"wrote {PI1_PI0_CSV}")

    pid_results = run_pid_convergence()
    _write_results(pid_results, PID_CSV)
    print(f"wrote {PID_CSV}")


if __name__ == "__main__":
    main()
