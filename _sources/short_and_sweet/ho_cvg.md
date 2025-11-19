# Convergence Rates


NG-BEM implements high order boundary element methods. The accuracy of the method depends on the order of the finite element spaces involved. On this page some numerical results of NG-BEM are presented: 

1. Analytically solvable problems with an unkown in one of the following trace spaces are considered 

$$ H^{\frac12}(\Gamma), \quad \boldsymbol H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma) \quad \textnormal{and} \quad H^{-\frac12}(\Gamma)\,.$$ 

The unkown trace is computed by NG-BEM for different orders of approximation $p$ on a sequence of subsequently refined meshes. The absolute $L_2$-error as function of degrees of freedom is analysed and convergence rates known from theory are exemplarily verified. 


2. Given a complete set of Cauchy data the solution of the pde is given by its boundary integral representation. Thus, besides traces the numerical solution of the pde depends on $p$. We consider an electrostatic potential in $H^{1}(\Omega)$ as solution of the Laplacian and analyse the decrease of error with respect to varying order $p$.


**Test 1: Laplace Dirichlet trace in $H^{\frac12}(\Gamma)$** 

Consider the harmonic function $u$,   

$$ u: \boldsymbol x \mapsto \dfrac{1}{\| \boldsymbol x - \boldsymbol x_s\|}, \quad \boldsymbol x_s = (1,1,1)^\intercal \in \Omega^c\,, $$ 

as unique solution of the Neumann boundary value problem

|  |  |
| -|--|
|$\begin{array}{rcl l} \Delta u &=& 0, \quad  & \mathrm{in} \, \Omega\,, \\ \gamma_1 u &=& u_1, \quad & \mathrm{on} \, \Gamma\,. \end{array}$    |  ![](resources/BEM_interior.png)  |

Let $\Omega$ be the unit ball and $\Gamma$ the unit sphere, respectively. Using the NG-BEM solver we compute the Dirichlet data for varying order $p$ on a sequence of meshes with decreasing mesh size, i.e., $h \sim\frac1N$. As the exact solution is kown, we can compute the absolute $L_2$-error of the numerical solution. This experiment is a proof of concept demonstrating examplarily the high order convergence rates that you can get with of a high order BEM solving for the Dirchlet trace:


![](resources/laplace_NtD_accuracy.png) 

Note that the Neumann problem is not uniquely solvable. They are the same up to a constant. Thus we compare here the surface gradients of the numerical and the exact solution.



**Test 2: Laplace Neumann trace in $H^{-\frac12}(\Gamma)$** 

To analyse convergence rates for traces in $H^{-\frac12}(\Gamma)$ consider now the Dirichlet boundary value problem with unkonwn Neumann data, i.e., 


|  |  |
| -|--|
|$\begin{array}{rcl l} \Delta u &=& 0, \quad  & \mathrm{in} \, \Omega\,, \\ \gamma_0 u &=& u_0, \quad & \mathrm{on} \, \Gamma\,. \end{array}$    |  ![](resources/BEM_interior.png)  |

Thus, again $u$ from Test 1 is the unique solution of this problem. Using the NG-BEM solver we compute the Neumann data for varying order $p$  on a sequence of meshes with decreasing mesh size, i.e., $h \sim\frac1N$. As the exact solution is kown, we compute the absolute $L_2$-error of the numerical solution. This experiment is a proof of concept demonstrating examplarily the high order convergence rates that you can get with of a high order BEM solving for the Dirchlet trace:



![](resources/laplace_DtN_accuracy.png) 


**Common notes for Test 1 and 2:** 

* Our manufactured solution $u$ does not depend on the mesh, i.e., it is a solution for all meshes. This is why the geometrical approximation order is $1$ for all test runs.
* Check out the repository's `convergence_timing` folder if you like to run this test.



**Test 3: Maxwell Neumann trace in $H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma)$** 

Besides electrostatics we can also solve boundary value problems from electromagnetics. The relevant trace spaces are $\boldsymbol H^{-\frac12}(\mathrm{curl}_\Gamma, \Gamma)$, and its dual, $\boldsymbol H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma)$. NG-BEM offers finite element spaces for both trace spaces. 

In order to verify convergence rates for high order bem, we consider the scattering of a plane wave $\boldsymbol E_{\mathrm{inc}}$ at a perfect electric conducting sphere. This scattering problem is analytically solvable by so-called Mie series. In order to keep the geometrical approximation error negligible, the geometrical approximation oder is set to $4$ for all tests runs. 

The scattered electric field $\boldsymbol E$ solves the following boundary value problem: 

|  |  |
| -|--|
|$\begin{array}{rcl l} \mathbf{\mathrm{curl}}\,\mathbf{\mathrm{curl}}\, \boldsymbol E - \kappa^2 \, \boldsymbol E&=& \boldsymbol 0, \quad  & \mathrm{in} \, \Omega^c\,, \\ \gamma_R \boldsymbol E &=& \boldsymbol E_{\mathrm{inc}}, \quad & \mathrm{on} \, \Gamma\,, \\ \textnormal{ + radiation  condition} & & &\|x\| \to \infty\,. \end{array}$    |  ![](resources/BEM_exterior.png)  |


Using the NG-BEM solver we compute the Neumann data for varying order $p$ on a sequence of meshes with decreasing mesh size i.e., $h \sim\frac1N$. The mesh is approximated with curvilinear elements of order $4$. As the exact solution for Mie-scattering kown, we compute the absolute $L_2$-error of the numerical solution.  This experiment is a proof of concept demonstrating examplarily the high order convergence rates that you can get with of a high order BEM solving for the surface current density:


![](resources/maxwell_mie_accuracy.png) 


**Notes:** 

* The Mie series solves the scattering at the sphere. A high order mesh is required to obtain high order convergence at least if the traces are approximated with $p>2$.   
* For approximation order $p=0$ the linear edge element functions are used, i.e., for a discretisation of trace space $H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma)$ those funcitons are called Rao-Wilton-Glisson functions (RWG). 
* Check out the repository's `convergence_timing` folder if you like to run this test.

