Laplace
-----------------------------

#### Layer Potentials

$$ \begin{array}{r rcl} \mathrm{SL}\left( j \right) (\boldsymbol x) &=& \displaystyle{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{1}{\| \boldsymbol x- \boldsymbol y\|} } \, j(\boldsymbol y)\, \mathrm{d}\sigma_y } \\ 
 \mathrm{DL}\left(m \right)(\boldsymbol x)  &=& \displaystyle{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{ \langle \boldsymbol n_y, \boldsymbol x- \boldsymbol y \rangle }{\| \boldsymbol x- \boldsymbol y\|^3} } \, m(\boldsymbol y)\, \mathrm{d}\sigma_y }
\end{array}$$

#### Laplace Dirichlet BVP

Let $u$ denote the electrostatic potential that arises under given Dirichlet boundary condition inside a source-free domain $\Omega \in \mathbb R^3$. Thus, $u$ solves the interior boundary value problem 

$$ \left\{ \begin{array}{rcl l} \Delta u &=& 0\,, \quad & \Omega \subset \mathbb R^3\,, \\ \gamma_0 u &=& m\,, \quad & \Gamma = \partial \Omega\,. \end{array} \right. $$
 

From here we can choose an direct or an indirect ansatz.  

**1. Direct Ansatz**  

$$ \begin{array}{r rcl }  
\textnormal{ansatz} & u &=& \mathrm{SL}(j) - \mathrm{DL}(m) \\
\textnormal{variational formulation } & \forall v\in H^{-\frac12}(\Gamma): \, \big\langle \gamma_0 \left(\mathrm{SL}(j)\right), v \big\rangle_{-\frac12} &=& \big\langle m, v\big\rangle_{-\frac12} + \big\langle \gamma_0 \left(\mathrm{DL}(m)\right), v \big\rangle_{-\frac12} \\ 
 \textnormal{discretisation} & \mathrm{V} \,\mathrm{j} &=& \left( \frac12 \,\mathrm{M} + \mathrm{K} \right) \, \mathrm{m} \\ 
\end{array}$$

**2. Indirect Ansatz** 

$$ \begin{array}{r rcl }  
\textnormal{ansatz} & u &=& \mathrm{SL}(j)  \\
\textnormal{variational formulation } & \forall v\in H^{-\frac12}(\Gamma): \, \big\langle \gamma_0 \left(\mathrm{SL}(j)\right), v \big\rangle_{-\frac12} &=& \big\langle m, v\big\rangle_{-\frac12} \\ 
 \textnormal{discretisation} & \mathrm{V} \, \mathrm{j} &=& \mathrm{M} \,\mathrm{m}  
\end{array} $$ 


#### Laplace Neumann BVP

Let $u$ denote the electrostatic potential that arises under given Neumann boundary condition inside a source-free domain $\Omega \in \mathbb R^3$. Thus, $u$ solves the boundary value problem

$$ \left\{ \begin{array}{rcl l} \Delta u &=& 0\,, \quad & \Omega \subset \mathbb R^3\,, \\ \gamma_1 u &=& j\,, \quad & \Gamma = \partial \Omega\,. \end{array} \right. $$


From here we can choose an direct or an indirect ansatz. 

**1. Direct Ansatz** 

$$ \begin{array}{r rcl }  
\textnormal{ansatz} & u &=& \mathrm{SL}(j) - \mathrm{DL}(m) \\
\textnormal{variational formulation }  & \forall v\in H^{\frac12}(\Gamma): \, \big\langle v, \gamma_1 \left(\mathrm{DL}(m)\right) \big\rangle_{-\frac12}  &=& \big\langle v, j\big\rangle_{-\frac12} - \big\langle v, \gamma_1 \left(\mathrm{SL}(j)\right) \big\rangle_{-\frac12} \\ 
\textnormal{discretisation} & \left( \mathrm{D} + \mathrm{S}\right) \mathrm{m} &=& \left( \frac12 \mathrm{M} - \mathrm{K}' \right) \, \mathrm{j}  
\end{array} $$ 

**2. Indirect Ansatz** 

$$ \begin{array}{r rcl }  
\textnormal{ansatz} & u &=&  \mathrm{DL}(m) \\
\textnormal{variational formulation } & \forall \quad v\in H^{\frac12}(\Gamma):\, \big\langle v, \gamma_1 \left(\mathrm{DL}(m)\right) \big\rangle_{-\frac12} &=& -\big\langle v, j\big\rangle_{-\frac12} \\ 
\textnormal{discretisation} & \left( \mathrm{D} + \mathrm S\right) \, \mathrm{m} &=&  -\mathrm{M}\,\mathrm{j}  
\end{array} $$ 


<!--
#### NG-BEM Potentials

Boundary integral equations are given based on a representation of the pde solution in terms of single and/or double layer potentials:

|  | trace space | type |  
|-|-|-|
| single layer potential | $H^{-\frac12}(\Gamma)$ | Neumann trace |
| double layer potential | $H^{\frac12}(\Gamma)$  | Dirichlet trace |

- NG-BEM implements the layper potentials based on conforming finite element spaces. 
- The finite element spaces are either natural traces of energy spaces:
  - The trace space $H^{\frac12}(\Gamma)$ is naturally given by $\gamma_0$`H1`.
  - The trace space $H^{-\frac12}(\Gamma)$ which is explicitely implemented as finite element (FE) space `SurfaceL2`. 

| Python interface | FE trace space  |   
|-|-|
|`LaplaceSL` |  `SurfaceL2` |
|`LaplaceDL` | $\gamma_0$ `H1` |


#### NG-BEM Potential Operators

The discretiszation of the boundary integral equations leads to the following layer potential operators:

|  | trial space | test space |  
|-|-|-|
| single layer potential operator | $H^{-\frac12}(\Gamma)$ | $H^{-\frac12}(\Gamma)$ |
| double layer potential operator | $H^{\frac12}(\Gamma)$  | $H^{-\frac12}(\Gamma)$ |
| hypersingular operator          | $H^{\frac12}(\Gamma)$  | $H^{\frac12}(\Gamma)$  |
| adjoint double layer potential operator | $H^{-\frac12}(\Gamma)$ | $H^{\frac12}(\Gamma)$  | 

- NG-BEM implements the layper potential operators based on conforming finite element spaces. 
- The finite element spaces are either natural traces of energy spaces:
  - The trace space $H^{\frac12}(\Gamma)$ is naturally given by $\gamma_0$`H1`.
  - The trace space $H^{-\frac12}(\Gamma)$ which is explicitely implemented as finite element (FE) space `SurfaceL2`. 

| Python interface | symbol |  FE trial space | FE test space |   
|-|:-:|-|-|
|`SingleLayerPotentialOperator` | $\mathrm V $ |  `SurfaceL2` | `SurfaceL2`|
|`DoubleLayerPotentialOperator` | $\mathrm K $ | $\gamma_0$ `H1` | `SurfaceL2` |
|`HypersingularOperator       ` | $\mathrm D$  | $\gamma_0$ `H1` | $\gamma_0$ `H1` |
|`DoubleLayerPotentialOperator` | $\mathrm K'$ | `SurfaceL2` | $\gamma_0$ `H1` |               
-->
