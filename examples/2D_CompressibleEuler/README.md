This is the  spectra of the 2D Linearized Euler equations
$$		
\begin{align}
0&=\partial_t \rho + \nabla \cdot \begin{pmatrix} \rho v_x \\ \rho v_y \end{pmatrix}\\
		\boldsymbol 0  &= \partial_t \begin{pmatrix} \rho v_x \\ \rho v_y \end{pmatrix} + \nabla \cdot \begin{pmatrix}
			\rho v_x^2 & \rho v_x v_y \\ \rho v_y v_x & \rho v_y^2 
		\end{pmatrix} + \nabla p \\
	  0 &= \partial_t E + \nabla \cdot \begin{pmatrix} (E + p) v_x \\ (E + p) v_y \end{pmatrix}
\end{align}$$
with total (internal + kinetic) energy
$$E = E(\rho, \boldsymbol v, p) = \rho \left(\frac{p}{\gamma - 1} + \frac{1}{2} \left(v_x^2 + v_y^2\right) \right).$$

We are interested in the advection of an isentropic vortex.

We discretize with the Discontinuous Galerkin Spectral Element Method (DGSEM) with
sixth order local polynomials and HLLC flux on $\Omega = [-10,10]^2$ with 64 elements in each direction.

Exemplary reference data:

First order accurate method:

`S = 16, t = 0.644182175756768579`

Second order accurate method:

`S = 16, t = 0.644076419407201705`
