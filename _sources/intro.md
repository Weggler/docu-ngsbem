# Welcome to NGSBEM

**NGSBEM** implements boundary integral operators on top of [NGSolve].  

It currently supports **single-layer, double-layer, and hypersingular operators** for:  

* Laplace equation  
* Helmholtz equation  
* Maxwell equation  

The software works with **high-order function spaces on curved surface meshes**.  
Numerical integration follows {cite}`SauterSchwab2009`, while **matrix compression and potential evaluation** are based on the **Fast Multipole Method** {cite}`RokhlinGreengard1985`.  

NGSBem achieves **high-order convergence rates** such as presented and discussed in {cite}`Weggler2011`.  

Install a [NGSolve release](https://github.com/NGSolve) and try notebooks from the [GitHub demos folder](https://github.com/Weggler/docu-ngsbem).  

---

##### Overview

The repository combines **practical demos** with a **theoretical introduction** to BEM:

#### Short and Sweet: Introduction to BEM and Software Capabilities

This section gives a concise overview of the **Boundary Element Method (BEM)** and the features of NGSBem.  
It highlights the **convergence rates** achievable for different problems.  

#### Demos: Step-by-Step from PDE to BEM Solution  

The demos show how to solve concrete problems using NGSBEM:  

* Explain the connection between **linear operators** and the **boundary value problem**.  
* Show step by step how the problem is translated into the boundary element method.  
* Introduce **software-specific operators** and key **implementation details**.  

These demos provide a hands-on view of the workflow, making the theory concrete.  

#### Background: Theoretical Foundations
This section introduces **energy spaces** and **trace spaces**.  

* Unlike FEM, BEM operates **on the boundary**, not the full PDE domain.  
* Function spaces are therefore discretized on the boundary.  
* All relevant **potentials** forming the starting point of BEM are presented.  

The presentation is theoretical and mainly follows {cite}`SauterSchwab2009` and {cite}`Weggler2011`.  

---

Overall, the repository offers a **hands-on introduction to NGSBEM**, linking practical demos to a concise summary of the **theoretical foundations**.  
Explore it on GitHub: [https://github.com/Weggler/docu-ngsbem](https://github.com/Weggler/docu-ngsbem).  

