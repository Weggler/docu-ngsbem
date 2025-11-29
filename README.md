# docu-ngsbem


Netgen/NGSolve is a high-performance finite element library designed for the numerical solution of partial differential equations. Developed in C++/C, it offers a streamlined Python interface for ease of use.
 
Since August 2025, NGSolve has incorporated boundary element method (BEM) functionality. The development of the NGSBem module was initiated in November 2023 by Lucy Weggler, who brings extensive experience with high-order boundary element methods.

This repository is intended as a supplement to the existing NGSolve tutorials and demos. 

- It presents a carefully curated collection of notebooks in the demos folder. The demos illustrate how NGSolve supports the implementation of boundary element methods. The examples demonstrate the breadth of current BEM capabilities, specifically showcasing the Python interface for solving boundary value problems related to the homogeneous Laplace, Helmholtz, and Maxwell equations as available in NGSolve.
- The Python scripts in the `convergence_timing` folder implement systematic convergence tests, specifically h-versions, for model problems with a known non-trivial analytical solution. Using the analytical solution, the exact error in the data can be computed, which the BEM provides for the respective trace.

For additional online documentation with demos and a touch of theory, visit: https://weggler.github.io/docu-ngsbem/intro.html

For NGSolve documentation, please see: https://docu.ngsolve.org/latest/


