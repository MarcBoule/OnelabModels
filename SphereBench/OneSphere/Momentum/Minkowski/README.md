# Minkowski Momentum Density

This folder contains simulation files for testing the Minkowski momentum density ($$\mathbf{D}\times\mathbf{B}$$). Two test cases from the parent folder are studied here: the calculation of the linear momentum in the "P_M" test case and the angular momentum in the "Rho_M" test case.

The open space integrals become (full domain $\Omega$):

$$
\mathbf{p}=\int_\Omega \mathbf{D} \times \mathbf{B}\,d\Omega,\qquad\quad
\mathbf{L}=\int_\Omega \mathbf{r} \times (\mathbf{D} \times \mathbf{B})\,d\Omega
$$

The finite space integrals become (sphere only $\Omega_s$):

$$
\mathbf{p}=
\varepsilon_0\mu_0\int_{\Omega_s} \mathbf{E} \times \mathbf{M}\,d\Omega_s
+\int_{\Omega_s} \mathbf{P} \times \mathbf{B}\,d\Omega_s
$$

$$
\mathbf{L}=
\int_{\Omega_s} \mathbf{r} \times (\rho\mathbf{A})\,d\Omega_s
+\int_{\Omega_s} \mathbf{r} \times (\mathbf{P}\times\mathbf{B})\,d\Omega_s
$$
