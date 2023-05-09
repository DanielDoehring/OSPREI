This is the  spectra of the 2D Linearized Euler equations
$$ 	\partial_t
		\begin{pmatrix}
			\rho' \\ u' \\ v' \\ p'
		\end{pmatrix}
		+
		\partial_x
		\begin{pmatrix}
			\bar{\rho} u' + \bar{u} \rho ' \\ \bar{u} u' + \frac{p'}{\bar{\rho}} \\ \bar{u} v' \\ \bar{u} p' + c^2 \bar{\rho} u'
		\end{pmatrix}
		+
		\partial_y
		\begin{pmatrix}
			\bar{\rho} v' + \bar{v} \rho ' \\ \bar{v} u' \\ \bar{v} v' + \frac{p'}{\bar{\rho}} \\ \bar{v} p' + c^2 \bar{\rho} v'
		\end{pmatrix}
		=
		\boldsymbol s $$
with source terms 
$$ 	\boldsymbol s(t,x,y) = 2\pi \begin{pmatrix}
			-\cos(2\pi t) \left[\cos(2\pi x) - \cos(2\pi y) \right] \\ \sin(2\pi t)\sin(2 \pi x) \\ \sin(2\pi t)\sin(2 \pi y) \\ -\cos(2\pi t) \left[\cos(2\pi x) - \cos(2\pi y) \right]
		\end{pmatrix}	
$$

corresponding to solution
$$ \begin{pmatrix}
			\rho'\\ u'\\ v' \\ p'
		\end{pmatrix}
		= \begin{pmatrix}
			-\cos(2\pi t) \left[\sin(2\pi x) - \sin(2\pi y) \right] \\ \sin(2\pi t)\cos(2 \pi x) \\ \sin(2\pi t)\cos(2 \pi y) \\ -\cos(2\pi t) \left[\sin(2\pi x) - \sin(2\pi y) \right]
		\end{pmatrix}.
$$

We discretize with the Discontinuous Galerkin Spectral Element Method (DGSEM) with
fourth order local polynomials and HLL flux on $\Omega = [0,1]^2$ with 64 elements in each direction.

Exemplary reference data:

First order accurate method:

`S = 16, t = 0.0457249149196286436`

Second order accurate method:

`S = 16, t = 0.0457002736631693534`
