Layer Potentials
-----------------------------

**Laplace**

$$ \begin{array}{r rcl} \mathrm{SL}\left( j \right) (\boldsymbol x) &=& \displaystyle{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{1}{\| \boldsymbol x- \boldsymbol y\|} } \, j(\boldsymbol y)\, \mathrm{d}\sigma_y }\,, \quad j \in H^{-\frac12}(\Gamma) \\ 
 \mathrm{DL}\left(m \right)(\boldsymbol x)  &=& \displaystyle{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{ \langle \boldsymbol n_y, \boldsymbol x-\boldsymbol y \rangle }{\| \boldsymbol x- \boldsymbol y\|^3} } \, m(\boldsymbol y)\, \mathrm{d}\sigma_y }\,,\quad m \in H^{\frac12}(\Gamma)
\end{array}$$

**Helmholtz** 

$$ \begin{array}{r rcl} \mathrm{SL}\left( j \right) (\boldsymbol x) &=& \displaystyle{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, \frac{e^{i\, \kappa \, \|\boldsymbol x-\boldsymbol y\|} }{\| \boldsymbol x- \boldsymbol y\|} } \, j(\boldsymbol y)\, \mathrm{d}\sigma_y }, \,, \quad j\in H^{-\frac12}(\Gamma) \\ 
 \mathrm{DL}\left(m \right)(\boldsymbol x)  &=& \displaystyle{ \int\limits_\Gamma \displaystyle{\frac{1}{4\,\pi}\, e^{i\,\kappa\,\|\boldsymbol x-\boldsymbol y|} \, \left( i\,\kappa \, \|\boldsymbol x-\boldsymbol y\| -1\right)\, \frac{ \langle \boldsymbol n_y, \boldsymbol x- \boldsymbol y \rangle }{\| \boldsymbol x- \boldsymbol y\|^3} } \, m(\boldsymbol y)\, \mathrm{d}\sigma_y }\,,\quad  m \in H^{\frac12}(\Gamma)
\end{array}$$


**Maxwell**

$$ \begin{array}{rcl} 
\mathrm{SL}\big(\boldsymbol j\big)(\boldsymbol x) &=& \kappa \, \displaystyle {\int\limits_\Gamma \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|\boldsymbol x- \boldsymbol y\|}}{\| \boldsymbol x-\boldsymbol y\|} \, \boldsymbol j(\boldsymbol y)\, \mathrm{d}\sigma_y + \frac{1}{\kappa} \nabla \int\limits_\Gamma \frac{1}{4\,\pi}\, \frac{e^{i\,\kappa\,\|\boldsymbol x- \boldsymbol y\|}}{\| \boldsymbol x - \boldsymbol y\|}  \, \mathrm{div}_\Gamma \boldsymbol j(\boldsymbol y)\, \mathrm{d}\sigma_y }\,,\quad \boldsymbol j\in H^{-\frac12}(\mathrm{div}_\Gamma, \Gamma) \\
\mathrm{DL}\big(\boldsymbol n \times \boldsymbol m\big)(\boldsymbol x)  &=& \nabla \times \displaystyle {\int\limits_\Gamma \displaystyle{ \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,\|\boldsymbol x-\boldsymbol y\|}}{\| \boldsymbol x- \boldsymbol y\|} } \, \boldsymbol n_y \times \boldsymbol{m}(\boldsymbol y)\, \mathrm{d}\sigma_y }\,, \quad  \boldsymbol m \in H^{-\frac12}(\mathrm{curl}_\Gamma, \Gamma)\end{array}
$$ 

**Trace Operators**

$$ \begin{array}{r rcl } \textnormal{Dirichlet trace:} \; & \gamma_0 u &=& u  \\[1ex] \textnormal{Neumann trace} \quad & \gamma_1 u &=& \langle \boldsymbol n,   \nabla\, u \rangle \,. \end{array} $$

- densities in $ H^{\frac12}\left( \Gamma\right) $ are weakly continous
- densities in $ H^{-\frac12}\left( \Gamma\right) $ are not continous

$$ \begin{array}{r rcl } \textnormal{Dirichlet trace} \quad & \gamma_R \boldsymbol u &=& \left( \boldsymbol n \times \boldsymbol u \right) \times \boldsymbol n \\[1ex] \textnormal{rotated Dirichlet trace} \quad & \gamma_D \boldsymbol u &=& \boldsymbol n \times \boldsymbol u \\ \textnormal{Neumann trace} \quad & \gamma_N \boldsymbol u &=& \dfrac{1}{\kappa} \boldsymbol n \times \mathbf{curl}\, \boldsymbol u\,,\quad \kappa = \omega\, \sqrt{\varepsilon\,\mu} \\ \textnormal{normal trace} \quad & \gamma_{\boldsymbol n} \boldsymbol u &=& \big\langle \boldsymbol n, \boldsymbol u \big\rangle\,. \end{array} $$

- tangential edge projection of densities in  $H^{-\frac12}\left( \mathrm{curl}_\Gamma,\Gamma\right)$ are weakly continous
- normal edge projection of densities in  $H^{-\frac12}\left( \mathrm{div}_\Gamma,\Gamma\right)$ are weakly continous
- densities in  $H^{-\frac12}\left( \Gamma\right)$ are not continous

**Energy Spaces and Trace Spaces**

$$
\begin{array}{rcccccc}
\textnormal{natural sequence:} &H^{\frac12}(\Gamma) & \xrightarrow{\nabla_{\Gamma}} & \boldsymbol{H}^{-\frac12}(\mathrm{curl}_{\Gamma},{\Gamma}) & \xrightarrow{\mathrm{curl}_{\Gamma}}& H^{-\frac12}({\Gamma})& \\[1ex]
&\gamma_0 \Big\uparrow && \gamma_R \Big\uparrow && \gamma_{\boldsymbol n} \Big\uparrow &\\[1ex]
\textnormal{energy spaces:} &H^1({\Omega}) & \xrightarrow{\nabla} & H(\mathbf{curl},{\Omega}) & \xrightarrow{\mathbf{curl}}& H(\mathrm{div},{\Omega}) & \xrightarrow{\mathrm{div}} \; L_2(\Omega) \\[1ex]
&\gamma_0 \Big\downarrow && \gamma_D \Big\downarrow && \gamma_{\boldsymbol n} \Big\downarrow &\\[1ex]
\textnormal{dual sequence:} & H^{\frac12}(\Gamma) & \xrightarrow{\mathbf{curl}_{\Gamma}} & \boldsymbol{H}^{-\frac12}(\mathrm{div}_{\Gamma},{\Gamma}) & \xrightarrow{\mathrm{div}_{\Gamma}}& H^{-\frac12}({\Gamma})& 
\end{array}
$$

Notes: 
- the natural sequence does not need any metric (geometric), it is **intrinsic** in a way
- the dual sequence needs metric (i.e., normal vector), it is **extrinsic** in a way

NGSolve implements conforming finite element spaces for all trace and energy spaces. The trace spaces from the natural sequence are natural in the sense that they are restrictions of the functions inside the domain on the boundary:  
- the natural trace space of $H^1(\Omega)$ is $H^{\frac12}(\Gamma)$ - in NGSolve it is the restriction of on the boundary. 
- the natural trace space of $H(\mathbf{curl},\Omega)$ is $H^{-\frac12}(\mathrm{curl}_\Gamma,\Gamma)$ - in NGSolve it is the vector component which is normal to the normal vector (projection on the tangent plane). 

The trace spaces for $H^{-\frac12}(\Gamma)$ and the $H^{-\frac12}(\mathrm{div}_\Gamma,\Gamma)$ are extra implementations in NGSolve. Note that NGSolve offers for some spaces a boolean flag `dual` that allows to switch to the extrinsic spaces (details see {cite}`Weggler2011`).

$$
\begin{array}{rcccccc}
\textnormal{natural sequence:} & \verb-H1.on.Bndry- & \xrightarrow{\nabla_{\Gamma}} &  \verb-Hcurl.on.Bndry- & \xrightarrow{\mathrm{curl}_{\Gamma}}& \verb-SurfaceL2-& \\[1ex]
&\gamma_0 \Big\uparrow && \gamma_R \Big\uparrow && \gamma_{\boldsymbol n} \Big\uparrow &\\[1ex]
\textnormal{energy spaces:} & \verb-H1-& \xrightarrow{\nabla} & \verb-Hcurl- & \xrightarrow{\mathbf{curl}}& \verb-Hdiv- & \xrightarrow{\mathrm{div}} \; \verb-L2- \\[1ex]
&&& && \\[1ex]
\textnormal{dual sequence:} & & & \verb-HdivSurface- & & \verb-SurfaceL2.dual-& 
\end{array}
$$
<!--
$$
\begin{array}{rcccccc}
\textnormal{natural sequence:} & \verb-H1.on.Bndry- & \xrightarrow{\nabla_{\Gamma}} &  \verb-Hcurl.on.Bndry- & \xrightarrow{\mathrm{curl}_{\Gamma}}& \verb-SurfaceL2-& \\[1ex]
&\gamma_0 \Big\uparrow && \gamma_R \Big\uparrow && \gamma_{\boldsymbol n} \Big\uparrow &\\[1ex]
\textnormal{energy spaces:} & \verb-H1-& \xrightarrow{\nabla} & \verb-Hcurl- & \xrightarrow{\mathbf{curl}}& \verb-Hdiv- & \xrightarrow{\mathrm{div}} \; \verb-L2- \\[1ex]
&\gamma_0 \Big\downarrow && \gamma_D \Big\downarrow && \gamma_{\boldsymbol n} \Big\downarrow &\\[1ex]
\textnormal{dual sequence:} &\verb-H1.on.Bndry.dual- & \xrightarrow{\mathbf{curl}_{\Gamma}} & \verb-HdivSurface- & \xrightarrow{\mathrm{div}_{\Gamma}}& \verb-SurfaceL2.dual-& 
\end{array}
$$
-->
