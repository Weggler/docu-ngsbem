#!/usr/bin/env python
# coding: utf-8

from ngsolve import *
from netgen.occ import *
import netgen.meshing as meshing
from ngsolve.krylovspace import CG, GMRes, MinRes
from ngsolve.bem import *
import time
import pandas as pd
from pathlib import Path

csv_path = Path("results.csv")
results = []


# Bottom sphere: Dirichlet boundary, top sphere: Neumann boundary
sp = Sphere((0, 0, 0), 1)
# mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.2)).Curve(4)

# h-refinement for fixed order
print("order ndofL2 ndofH1 err_dir_semi time type")
for order in range(1, 5):
    for i in range(1, 7):
        # mesh = Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.75 / i, perfstepsend=meshing.MeshingStep.MESHSURFACE)) # .Curve(order)
        mesh = Mesh(OCCGeometry(sp).GenerateMesh(maxh=0.75 / i)).Curve(4)

        fesL2 = SurfaceL2(mesh, order=order, dual_mapping=False)
        u, v = fesL2.TnT()
        fesH1 = H1(mesh, order=order, definedon=mesh.Boundaries(".*"))
        uH1, vH1 = fesH1.TnT()

        uexa = 1 / sqrt((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2)

        graduexa = CF((uexa.Diff(x), uexa.Diff(y), uexa.Diff(z)))
        n = specialcf.normal(3)
        u1 = GridFunction(fesL2)
        u1.Interpolate(graduexa * n, definedon=mesh.Boundaries(".*"))

        intorder = 2 * order + 2  # +12
        vH1m1 = LinearForm(vH1 * 1 * ds(bonus_intorder=2 * intorder)).Assemble()
        S = (BaseMatrix(Matrix(vH1m1.vec.Reshape(1)))) @ (
            BaseMatrix(Matrix(vH1m1.vec.Reshape(fesH1.ndof)))
        )

        u0 = GridFunction(fesH1)
        start = time.time()
        with TaskManager():
            pre = (
                BilinearForm(uH1 * vH1 * ds(bonus_intorder=intorder))
                .Assemble()
                .mat.Inverse(freedofs=fesH1.FreeDofs())
            )
            D = LaplaceSL(curl(uH1) * ds(bonus_intorder=intorder)) * curl(vH1) * ds
            M = BilinearForm(v.Trace() * uH1 * ds(bonus_intorder=intorder)).Assemble()
            Kt = LaplaceDL(u * ds(bonus_intorder=intorder)) * vH1 * ds
            rhs = ((0.5 * M.mat.T - Kt.mat) * u1.vec).Evaluate()

            CG(
                mat=D.mat + S,
                pre=pre,
                rhs=rhs,
                sol=u0.vec,
                tol=1e-8,
                maxsteps=200,
                initialize=False,
                printrates=False,
            )

        end = time.time()
        elapsed = end - start

        u0exa = GridFunction(fesH1)
        u0exa.Interpolate(uexa, definedon=mesh.Boundaries(".*"))

        # Compute the L2-seminorm error
        errd = sqrt(
            Integrate((grad(u0exa) - grad(u0)) ** 2, mesh.Boundaries(".*"), BND)
        )

        print(order, fesL2.ndof, fesH1.ndof, errd, elapsed, "dirichlet_semi")

        results.append(
            {
                "order": order,
                "ndof": fesL2.ndof,
                "err": float(errd),
                "time": elapsed,
                "type": "dirichlet_semi",
            }
        )

df = pd.DataFrame(results)
if csv_path.exists():
    df.to_csv(csv_path, mode="a", header=False, index=False)
else:
    df.to_csv(csv_path, mode="w", header=True, index=False)
