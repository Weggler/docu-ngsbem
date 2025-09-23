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
sp = Sphere( (0,0,0), 1)
# mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.2)).Curve(4)

# h-refinement for fixed order
print("order ndofL2 ndofH1 err_neu time type")
for order in range(1, 5):
    for i in range(1, 7):
        # mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.75 / i, perfstepsend=meshing.MeshingStep.MESHSURFACE)).Curve(order)
        mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.75 / i))#.Curve(4)

        fesL2 = SurfaceL2(mesh, order=order-1, dual_mapping=True)
        u,v = fesL2.TnT()
        fesH1 = H1(mesh, order=order)
        uH1,vH1 = fesH1.TnT()
        
        uexa = 1/ sqrt( (x-1)**2 + (y-1)**2 + (z-1)**2 )
        u0 = GridFunction(fesH1)
        u0.Interpolate (uexa)
        u1 = GridFunction(fesL2)
        
        intorder = 2 * order + 2  # +12
        start = time.time()
        with TaskManager(): 
            pre = BilinearForm(u*v*ds, diagonal=True).Assemble().mat.Inverse()
            V = LaplaceSL(u*ds(bonus_intorder=intorder))*v*ds
            M = BilinearForm( uH1 * v * ds(bonus_intorder=intorder)).Assemble()
            K = LaplaceDL(uH1*ds(bonus_intorder=intorder))*v*ds
            
            rhs = ( (0.5 * M.mat + K.mat)*u0.vec).Evaluate()
            CG(mat = V.mat, pre=pre, rhs = rhs, sol=u1.vec, tol=1e-8, maxsteps=50, initialize=False, printrates=False)
            
        end = time.time()
        elapsed = end - start
        
        graduexa = CF( (uexa.Diff(x), uexa.Diff(y), uexa.Diff(z)) )
        n = specialcf.normal(3)
        u1exa = graduexa*n 

        # Compute the L2-norm error 
        errn = sqrt(Integrate((u1exa - u1)**2, mesh.Boundaries(".*"), BND))
        
        print(order, fesL2.ndof, fesH1.ndof, errn, elapsed, "neumann")
        
        results.append({
            "order": order,
            "ndof": fesL2.ndof,
            "err": float(errn),
            "time": elapsed,
            "type": "neumann"
        })

df = pd.DataFrame(results)
if csv_path.exists():
    df.to_csv(csv_path, mode="a", header=False, index=False)
else:
    df.to_csv(csv_path, mode="w", header=True, index=False)
