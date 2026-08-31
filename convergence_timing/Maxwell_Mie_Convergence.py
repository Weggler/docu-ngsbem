#!/usr/bin/env python

from netgen.occ import *
import netgen.meshing as meshing
from ngsolve import *
from ngsolve.bem import *
from ngsolve.krylovspace import GMRes
from ngsolve.fem import CompilePythonModule
from pathlib import Path
import sys
import time

from convergence_timing.common import append_results, max_mesh_size

csv_path = Path("bem_results.csv")
miecurrent = None

def get_mie_current():
    global miecurrent

    if miecurrent is None:
        source = "mie_ngs.cpp" if sys.platform == "darwin" else "mie.cpp"
        source_path = Path(__file__).with_name(source)

        if not source_path.exists():
            raise FileNotFoundError(
                f"Mie implementation not found for platform {sys.platform!r}: "
                f"{source_path}"
            )

        txt = source_path.read_text()

        mie = CompilePythonModule(
            txt,
            init_function_name="Mie",
            add_header=False,
        )

        miecurrent = mie.MieCurrent()

    return miecurrent

# Scattering on a sphere
sp = Sphere((0, 0, 0), 0.25)
kappa = 5.0
E_inc = CF((1, 0, 0)) * exp(1j * kappa * z)


def run_maxwell_mie(order, n):
    reference_current = get_mie_current()
    requested_maxh = 0.2 / n
    mesh = Mesh(
        OCCGeometry(sp).GenerateMesh(
            maxh=requested_maxh, perfstepsend=meshing.MeshingStep.MESHSURFACE
        )
    ).Curve(4)
    h = max_mesh_size(mesh)
    fesHDiv = HDivSurface(mesh, order=order, complex=True)
    uHDiv, vHDiv = fesHDiv.TnT()

    rhs = LinearForm(-E_inc * vHDiv.Trace() * ds(bonus_intorder=10)).Assemble()

    j = GridFunction(fesHDiv)
    start = time.time()
    intorder = order + 1
    with TaskManager():
        pre = (
            BilinearForm(
                uHDiv.Trace() * vHDiv.Trace() * ds(bonus_intorder=intorder)
            )
            .Assemble()
            .mat.Inverse(freedofs=fesHDiv.FreeDofs())
        )
        # V = MaxwellSingleLayerPotentialOperator(fesHDiv, kappa, intorder=intorder)
        # GMRes(A=V.mat, pre=pre, b=rhs.vec, x=j.vec, tol=1e-11, maxsteps=2000, printrates=False)
        V1 = (
            HelmholtzSL(uHDiv.Trace() * ds(bonus_intorder=intorder), kappa)
            * vHDiv.Trace()
            * ds(bonus_intorder=intorder)
        )
        V2 = (
            HelmholtzSL(div(uHDiv.Trace()) * ds(bonus_intorder=intorder), kappa)
            * div(vHDiv.Trace())
            * ds(bonus_intorder=intorder)
        )
        V = kappa * V1.mat - (1 / kappa) * V2.mat
        GMRes(
            A=V,
            pre=pre,
            b=rhs.vec,
            x=j.vec,
            tol=1e-11,
            maxsteps=2000,
            printrates=False,
        )

    end = time.time()
    elapsed = end - start

    j.vec[:] *= kappa

    error = sqrt(Integrate(Norm(j - reference_current) ** 2, mesh, BND))
    print(order, fesHDiv.ndof, error)

    return {
        "order": order,
        "h": h,
        "ndof": fesHDiv.ndof,
        "err": float(error),
        "time": elapsed,
        "type": "mie",
    }


def main():
    results = []

    print("order ndof error")
    for order in range(1, 5):
        for n in range(1, 7):
            results.append(run_maxwell_mie(order, n))

    append_results(results, csv_path)


if __name__ == "__main__":
    main()
