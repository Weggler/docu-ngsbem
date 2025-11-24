Boundary Element Method
------------------------

**What is a standard BEM problem?**

One of the following pdes with given boundary condition:

- Laplace equation 
- Helmholtz equation 
- Maxwells equations 
- Lamé equations 
- Stokes equations 

**How to derive a BEM from a standard BEM problem?**

1. choose an ansatz for the solution of the pde in terms of so called layer potentials
2. derive a boundary integral equation for unknown density
3. discretize the resulting variational formulation with finite element spaces on the boundary 
4. solve the system of linear equations and get the best approximation of the unknown density
4. evaluate the solution with the ansatz from 1. wherever you want inside the pde domain

**Why is BEM beneficial?** 

- problem dimension is reduced to the boundary of the pde domain, thus reduced by one
- exterior problems are not an issue
- the solution is very accurate 

**Why is BEM not everybodies darling?** 

- only linear, isotropic material 
- source terms cause Newton potentials 
- singular integral kernels
- dense matrices

**NGSBem Keyfacts:** 

- kernel-driven generic implementation of the layer potential operator 
- fast and accurate assembly of system matrices using MLFMM
- fast and accurate assembly of potentials using MLFMM with nearfield considerations
- compatible with NGSolve 
- user-friendly Python interface with declarative programming of variational formulations
- potentials for all standard problems
- tested and documented 

