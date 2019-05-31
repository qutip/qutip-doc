.. QuTiP
   Copyright (C) 2011-2012, Paul D. Nation & Robert J. Johansson

.. _stochastic:

*******************************************
Stochastic Solver
*******************************************

.. _stochastic-intro:


Introduction
=============





ssesolve
=========
.. Stochastic Schrodinger equation

.. math::
	:label: jump_ssesolve

    e =\left<\psi(t)|C + C^{+}|\psi(t)\right>\\
	\delta \psi(t) = - i H \psi(t) \delta t
	                 - \left(  \frac{C^{+} C}{2} -\frac{e}{2}C + \frac{e^2}{8} \right) \psi  \delta t
	                 + \left(  C -\frac{e}{2} \right) \psi  \delta \omega



smesolve
=========
.. Stochastic Master equation

.. math::
	:label: master_equation

	L(\rho(t))= - i[H(t),\rho(t)]
	            + \sum_n \frac{1}{2} \left[2 S \rho(t) S^{+} - \rho(t) S^{+} S - S^{+} S \rho(t)\right]
	            + \sum_n \frac{1}{2} \left[2 C \rho(t) C^{+} - \rho(t) C^{+} C - C^{+} C \rho(t)\right]

.. math::
	:label: jump_smesolve

	\delta \rho(t) = L(\rho(t)) \delta t
	                 + \left(  C \rho(t) + \rho(t) C^{+} - tr\left(C \times \rho + \rho \times C^{+} \right)\rho(t) \right)  \delta \omega



===============================================================================
