This is the  spectra of Burgers equation 
$$ u_t + \frac{1}{2} (u^2)_x = 0 $$
when discretized with the Discontinuous Galerkin Spectral Element Method (DGSEM).
Third order local polynomials and the Godunov flux are employed.
The jacobian is evaluated at $u_0(x) = 2 + \sin(2 \pi x)$.

Exemplary reference data:

First order accurate method:

`S = 16, t = 0.00294916495797224349`

Second order accurate method:

`S = 16, t = 0.00291312662768177694`
