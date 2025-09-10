#!/usr/bin/env python
# coding: utf-8

from ngsolve import *
from netgen.occ import *
import netgen.meshing as meshing
from ngsolve.krylovspace import CG, GMRes
from ngsolve.bem import *

# Bottom sphere: Dirichlet boundary, top sphere: Neumann boundary
topsphere = Sphere((0,0,0), 1) * Box((-1,-1,0),(1,1,1))
botsphere = Sphere((0,0,0), 1) - Box((-1,-1,0),(1,1,1))
topsphere.faces.name = "neumann"
botsphere.faces.name = "dirichlet"
shape = Fuse([topsphere,botsphere])
screen = WorkPlane(Axes((0,0,0), Z, X)).RectangleC(0.5,0.5).Face()
screen.faces.name="screen"
mesh_screen = Mesh(OCCGeometry(screen).GenerateMesh(maxh=0.3)).Curve(1)

print("order ndofL2 ndofH1 err_neu err_dir err_scr")
for order in range(1, 5):
    for i in range(1, 7):
        mesh = Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.75 / i, perfstepsend=meshing.MeshingStep.MESHSURFACE)) # .Curve(order)

        fesH1 = H1(mesh, order=order, dirichlet="dirichlet", definedon=mesh.Boundaries(".*"))
        fesL2 = SurfaceL2(mesh, order=order-1, dirichlet="neumann")
        fes = fesH1 * fesL2
        u,dudn = fes.TrialFunction()
        v,dvdn = fes.TestFunction()

        uexa = CF(1. / sqrt((x-1)**2 + (y-1)**2 + (z-1)**2))
        n = specialcf.normal(3)
        gradn_uexa = CF((uexa.Diff(x), uexa.Diff(y), uexa.Diff(z))) * n

        udn = GridFunction(fes)
        # Dirichlet data on Dirichlet boundary
        udn.components[0].Interpolate(uexa, definedon=mesh.Boundaries("dirichlet"))
        # Neumann data on Neumann boundary
        udn.components[1].Interpolate(gradn_uexa, definedon=mesh.Boundaries("neumann"))

        intorder = 2 * order + 2
        with TaskManager():
            V = LaplaceSL(dudn*ds(bonus_intorder=intorder)) * dvdn *ds
            K = LaplaceDL(u*ds(bonus_intorder=intorder)) * dvdn *ds
            W = LaplaceSL(curl(u)*ds(bonus_intorder=intorder)) * curl(v)*ds
            M = BilinearForm(u.Trace() * dvdn.Trace() * ds(bonus_intorder=intorder)).Assemble()

            lhs = V.mat - K.mat + K.mat.T + W.mat 
            rhs  = ((0.5 * ( M.mat + M.mat.T) + ( K.mat - K.mat.T) - V.mat - W.mat) * udn.vec).Evaluate()
            pre = BilinearForm( (u * v + dudn * dvdn) * ds(bonus_intorder=ntorder) ).Assemble().mat.Inverse(freedofs=fes.FreeDofs())

            sol = GMRes(A=lhs, b=rhs, pre=pre, maxsteps=400, printrates=False, tol=1e-12)

        gfu = GridFunction(fes)
        gfu.vec[:] = sol

        # Compute the L2-error in the traces on the complete boundary
        errd = sqrt(Integrate((uexa - gfu.components[0] - udn.components[0])**2, mesh.Boundaries(".*"), BND))
        errn = sqrt(Integrate((gradn_uexa - gfu.components[1] - udn.components[1])**2, mesh.Boundaries(".*"), BND))

        # Post-processing on screen
        fes_screen = H1(mesh_screen, order=13)
        gf_screen = GridFunction(fes_screen)
        uH1, vH1 = fesH1.TnT()
        uL2, vL2 = fesL2.TnT()
        SLPotential = LaplaceSL(uL2*ds(bonus_intorder=intorder))
        DLPotential = LaplaceDL(uH1*ds(bonus_intorder=intorder))
        repformula = SLPotential(gfu.components[1]) - DLPotential(gfu.components[0]) + SLPotential(udn.components[1]) - DLPotential(udn.components[0])
        with TaskManager():
            gf_screen.Set(repformula, definedon=mesh_screen.Boundaries("screen"), dual=False)

        errs = sqrt(Integrate((gf_screen - uexa)**2, mesh_screen.Boundaries(".*"), BND))
        # compute error on screen, without interpolation:
        errsint = sqrt(Integrate((repformula - uexa)**2, mesh_screen.Boundaries(".*"), BND))
        errsintinterpol = sqrt(Integrate((repformula - uexa)**2, mesh_screen.Boundaries(".*"), BND))

        print(order, fesL2.ndof, fesH1.ndof, errn, errd, errs, errsint, errsintinterpol)
