Helmholtz
-------------------------------

<!--
**Notations of trace operators:**

$$ \begin{array}{r rcl } \textnormal{Dirichlet trace} \quad & \gamma_0 u &=& u  \\[1ex] \textnormal{Neumann trace} \quad & \gamma_1 u &=& \langle \boldsymbol n,   \nabla\, u \rangle \,. \end{array} $$

**Properties of trace spaces:**

$$ \begin{array}{r rcl l} \textnormal{Dirichlet trace} \quad & \gamma_0 u &\in& H^{\frac12}\left( \Gamma\right) \quad &\textnormal{weakly continuous}\\  \textnormal{Neumann trace} \quad & \gamma_1 u &\in& H^{-\frac12}\left( \Gamma\right) \quad & \textnormal{not continuous}\,. \end{array} $$
-->
#### Layer Potentials

$$ \begin{array}{rcl} \mathrm{SL}\left( j \right) (x) &=& \displaystyle{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{e^{i\, \kappa \, |x-y|} }{\| x-y\|} } \, j(y)\, \mathrm{d}\sigma_y } \\ 
 \mathrm{DL}\left(m \right)(x)  &=& \displaystyle{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \displaystyle{ n_y \cdot \nabla_y }\displaystyle{ \frac{ e^{i\,\kappa\,\|x-y\|}}{\| x-y\|}} } \, m(y)\, \mathrm{d}\sigma_y }
\end{array}$$


#### Helmholtz Dirichlet BVP

Let $u$ denote the acoustic potential which is caused by Dirichlet boundary condition on $\Gamma$ and which propagates in $\Omega^c \in \mathbb R^3$. Thus, $u$ solves the exterior boundary value problem 

$$ \left\{ \begin{array}{rcl l} \Delta u + \kappa^2 u &=& 0\,, \quad & \Omega^c \subset \mathbb R^3\,, \\ \gamma_0 u &=& m\,, \quad & \Gamma = \partial \Omega\,. \end{array} \right. $$
 

To stabilize interior eigenvalues, one consideres a combined field integral equation. The combinded field integral equation combines single and double layer integral operators, one option is the **Brakhage-Werner** formulation: 

$$ \begin{array}{ll rcl }  
&\textnormal{ansatz} & u(x) &=& i \, \kappa \,  \mathrm{ SL}(j) - \mathrm{DL}(j) \\ 
&\textnormal{variational formulation }  &  \forall v\in H^{-\frac12}(\Gamma): \; \left\langle \gamma_0 u , v \right\rangle_{-\frac12} &=& i \, \kappa \, \left\langle \gamma_0 \left(\mathrm{SL}(j)\right), v \right\rangle_{-\frac12} - \left\langle \gamma_0 \left(\mathrm{DL}(j)\right), v\right\rangle_{-\frac12} \\ 
& \textnormal{discretisation} & \left( \frac12 \, \mathrm M + i \, \kappa \,\mathrm{V} + \mathrm K\right) \, \mathrm{j} &=& \mathrm{M} \, \mathrm{m} 
 \end{array} $$ 

<!--
#### NG-BEM Python Functions 

| Symbol | Operator | trial space | test space | NG-BEM | trial NG-Solve | test NG-Solve |   
|-|-|-|-|-|-|-|
| $\mathrm V $ | single layer | $H^{-\frac12}(\Gamma)$ | $H^{-\frac12}(\Gamma)$         | HelmholtzSingleLayerPotentialOperator | SurfaceL2 | SurfaceL2|
| $\mathrm K $ | double layer | $H^{\frac12}(\Gamma)$ | $H^{-\frac12}(\Gamma)$          | HelmholtzDoubleLayerPotentialOperator | H1 | SurfaceL2 |
| $\mathrm D$ | hypersingular  | $H^{\frac12}(\Gamma)$ | $H^{\frac12}(\Gamma)$          | HelmholtzHypersingularOperator | H1 | H1 |
| $\mathrm K' $ | adjoint double layer | $H^{-\frac12}(\Gamma)$ | $H^{\frac12}(\Gamma)$ | HelmholtzDoubleLayerPotentialOperator | SurfaceL2 | H1 |               
-->
