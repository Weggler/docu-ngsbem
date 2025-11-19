# Welcome to NGSBEM

NGSBem implements boundary integral operators on top of NGSolve.

Currently it supports single-layer and double-layer operators and hypersingular operators for the

* Laplace equation
* Helmholtz equation
* Maxwell equation

NGSBem supports high order function spaces on curved surface meshes.
It uses numerical integration following [Sauter-Schwab: Boundary Element Methods, 2009], and matrix compression and potential evalution based on the Fast Multipole Method  [Rokhlin, Greengard, Rapid Solution of Integral Equations of Classic Potential Theory, 1985].


Installation: install a recent NGSolve (later than November 19, 2025)

Try notebooks from the demos folder.
You can read (only read) them directly on github: https://github.com/Weggler/docu-ngsbem


```{tableofcontents}
```
