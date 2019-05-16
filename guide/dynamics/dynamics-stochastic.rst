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



photocurent stochastic equation
================================
.. photocurent Schrodinger equation


.. math::
	:label: jump_pssesolve

	\delta \omega = poisson \left(\left<\psi(t)|C^{+}C|\psi(t)\right> \right)\\

.. math::
	:label: jump_pssesolve2

	\delta \psi(t) = - i H \psi(t) \delta t
					 -\frac{C^{+} C}{2}  \psi(t) \delta t
					 + \frac{ \left| C \psi  \right| ^2}{2} \delta t
	                 + \delta \omega \left( \frac{C \psi}{\left| C \psi  \right|} - \psi \right)\\


photocurent Master equation
================================
.. photocurent Master equation

.. math::
	:label: master_equation

	L(\rho(t))= - i[H(t),\rho(t)]
	            + \sum_n \frac{1}{2} \left[2 S \rho(t) S^{+} - \rho(t) S^{+} S - S^{+} S \rho(t)\right]

.. math::
	:label: jump_psmesolve

	\delta \omega = poisson \left(tr \left( C \rho C^{+} \right) \right)\\

	\delta \rho(t) = L(\rho) \delta t  + \left( tr \left(C^{+}C  \rho C^{+}C \right) - C^{+}C  \rho C^{+}C \right) \frac{\delta t}{2}
	                 +\delta \omega \left( \frac{C \rho C^{+}}{tr \left( C \rho C^{+} \right)} - \rho \right)







===========================================================================================================================================
