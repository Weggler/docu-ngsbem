#!/usr/bin/env python

from netgen.occ import *
import netgen.meshing as meshing
from ngsolve import *
from ngsolve.bem import *
from ngsolve.krylovspace import GMRes
from ngsolve.fem import CompilePythonModule
from pathlib import Path
import time
import pandas as pd

csv_path = Path("results.csv")
results = []

# Get reference mie series current from cpp
txt = Path("mie_ngs.cpp").read_text()
mie = CompilePythonModule(txt, init_function_name="Mie", add_header=False)
miecurrent = mie.MieCurrent()

# Scattering on a sphere
sp = Sphere((0, 0, 0), 0.25)
kappa = 5.0
E_inc = CF((1, 0, 0)) * exp(1j * kappa * z)

print("order ndof error")
for order in range(1, 5):
    for n in range(1, 7):
        h = 0.2 / n
        mesh = Mesh(
            OCCGeometry(sp).GenerateMesh(
                maxh=h, perfstepsend=meshing.MeshingStep.MESHSURFACE
            )
        ).Curve(4)
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

        error = sqrt(Integrate(Norm(j - miecurrent) ** 2, mesh, BND))
        print(order, fesHDiv.ndof, error)

        results.append(
            {
                "order": order,
                "ndof": fesHDiv.ndof,
                "err": float(error),
                "time": elapsed,
                "type": "mie",
            }
        )

df = pd.DataFrame(results)
if csv_path.exists():
    df.to_csv(csv_path, mode="a", header=False, index=False)
else:
    df.to_csv(csv_path, mode="w", header=True, index=False)
