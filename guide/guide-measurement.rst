.. QuTiP
   Copyright (C) 2011-2012, Paul D. Nation & Robert J. Johansson

.. _measurement:

******************************
Measurement of Quantum Objects
******************************

.. ipython::
   :suppress:

   In [1]: from qutip import basis, sigmax, sigmaz

   In [2]: from qutip.measurement import measure, measurement_statistics


.. _measurement-intro:

Introduction
============

Measurement is a fundamental part of the standard formulation of quantum
mechanics and is the process by which classical readings are obtained from
a quantum object. Although the intrepretation of the procedure is at times
contentious, the procedure itself is mathematically straightforward and is
described in many good introductory texts.

Here we will show how to perform simple measurement operations on QuTiP
objects.

.. _measurement-basic:

Performing a basic measurement
------------------------------

First we need to select some states to measure. For now, let us create an *up*
state and a *down* state:

.. ipython::

   In [1]: up = basis(2, 0)

   In [2]: down = basis(2, 1)

which represent spin-1/2 particles with their spin pointing either up or down
along the z-axis.

We choose what to measure by selecting a *measurement operator*. For example,
we could select :func:`~qutip.sigmaz` which measures the z-component of the
spin of a spin-1/2 particle, or :func:`~qutip.sigmax` which measures the
x-component:

.. ipython::

   In [3]: spin_z = sigmaz()

   In [4]: spin_x = sigmax()

How do we know what these operators measure? The answer lies in the measurement
procedure itself:

* A quantum measurement tranforms the state being measure by projecting it into
  one of the eigenvectors of the measurement operator.

* Which eigenvector to project onto is chosen probabilitically according to the
  square of the amplitude of the state in the direction of the eigenvector.

* The value returned by the measurement is the eigenvalue corresponding to the
  chosen eigenvector.

.. note::

   How to intrepret this "random chosing" is the famous
   "quantum measurement problem".

The eigenvectors of `spin_z` are the states with their spin pointing either up
or down, so it measures the component of the spin along the z-axis.

The eigenvectors of `spin_x` are the states with their spin pointing either
left or right, so it measures the component of the spin along the x-axis.

When we measure our `up` and `down` states using the operator `spin_z`, we
always obtain:

.. ipython::

   In [5]: measure(spin_z, up) == (1.0, up)

   In [6]: measure(spin_z, down) == (-1.0, down)

because `up` is the eigenvector of `spin_z` with eigenvalue `1.0` and `down`
is the eigenvector with eigenvalue `-1.0`.

Neither eigenvector has any component in the direction of the other (they are
orthogonal), so `measure(spin_z, up)` returns the state `up` 100% percent of the
time and `measure(spin_z, down)` returns the state `down` 100% of the time.

Note how :func:`~qutip.measurement.measure` returns a pair of values. The
first is the measured value, i.e. an eigenvalue of the operator (e.g. `1.0`),
and the second is the state of the quantum system after the measurement,
i.e. an eigenvector of the operator (e.g. `up`).

Now let us consider what happens if we measure the x-component of the spin
of `up`:

.. ipython::

   In [7]: measure(spin_x, up)

The `up` state is not an eigenvector of `spin_x`. `spin_x` has two eigenvectors
which we will call `left` and `right`. The `up` state has equal components in
the direction of these two vectors, so measurement will select each of them
50% of the time.

When `left` is chosen, the result of the measurement will be `(1.0, left)`.

When `right` is chosen, the result of measurement with be `(-1.0, right)`.

Now you know how to measure quantum states in QuTip!

The `measure` function can perform measurements on density matrices too. You
can read about these and other details at :func:`~qutip.measurement.measure`.

.. _measurement-statistics:

Obtaining measurement statistics
--------------------------------

XXX: TODO
