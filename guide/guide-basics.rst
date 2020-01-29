.. QuTiP
   Copyright (C) 2011-2012, Paul D. Nation & Robert J. Johansson

.. _basics:

************************************
Basic Operations on Quantum Objects
************************************

.. _basics-first:

First things first
==================

.. warning:: Do not run QuTiP from the installation directory.

To load the qutip modules, we must first call the import statement:

.. ipython::

    In [1]: from qutip import *


that will load all of the user available functions. Often, we also need to import the NumPy and Matplotlib libraries with:

.. ipython::

   In [1]: import numpy as np

   In [2]: import matplotlib.pyplot as plt

Note that, in the rest of the documentation, functions are written using `qutip.module.function()` notation which links to the corresponding function in the QuTiP API: :ref:`functions`. However, in calling `import *`, we have already loaded all of the QuTiP modules. Therefore, we will only need the function name and not the complete path when calling the function from the interpreter prompt, Python script, or Jupyter notebook.

.. _basics-qobj:

The quantum object class
========================

.. _basics-qobj-intro:

Introduction
---------------

The key difference between classical and quantum mechanics lies in the use of operators instead of numbers as variables. Moreover, we need to specify state vectors and their properties. Therefore, in computing the dynamics of quantum systems we need a data structure that is capable of encapsulating the properties of a quantum operator and ket/bra vectors. The quantum object class, :func:`qutip.Qobj`, accomplishes this using matrix representation.

To begin, let us create a blank ``Qobj``:

.. ipython::

    In [1]: Qobj()

where we see the blank ``Qobj`` object with dimensions, shape, and data. Here the data corresponds to a 1x1-dimensional matrix consisting of a single zero entry.

.. Hint:: By convention, Class objects in Python such as ``Qobj()`` differ from functions in the use of a beginning capital letter.

We can create a ``Qobj`` with a user defined data set by passing a list or array of data into the ``Qobj``:

.. ipython::

    In [1]: Qobj([[1],[2],[3],[4],[5]])

    In [2]: x = np.array([[1, 2, 3, 4, 5]])

    In [3]: Qobj(x)

    In [4]: r = np.random.rand(4, 4)

    In [5]: Qobj(r)


Notice how both the dims and shape change according to the input data.  Although dims and shape appear to have the same function, the difference will become quite clear in the section on :ref:`tensor products and partial traces <tensor>`.

.. note:: If you are running QuTiP from a python script you must use the :func:`print` function to view the Qobj attributes.

.. _basics-qobj-states:

States and operators
---------------------

Manually specifying the data for each quantum object is inefficient. Even more so when most objects correspond to commonly used types such as the ladder operators of a harmonic oscillator, the Pauli spin operators for a two-level system, or state vectors such as Fock states. Therefore, QuTiP includes predefined objects for a variety of states:

.. cssclass:: table-striped

+--------------------------+----------------------------------+----------------------------------------+
| States                   | Command (# means optional)       | Inputs                                 |
+==========================+==================================+========================================+
| Fock state ket vector    | ``basis(N,#m)``/``fock(N,#m)``   | N = number of levels in Hilbert space, |
|                          |                                  | m = level containing excitation        |
|                          |                                  | (0 if no m given)                      |
+--------------------------+----------------------------------+----------------------------------------+
| Fock density matrix      | ``fock_dm(N,#p)``                | same as basis(N,m) / fock(N,m)         |
| (outer product of basis) |                                  |                                        |
+--------------------------+----------------------------------+----------------------------------------+
| Coherent state           | ``coherent(N,alpha)``            | alpha = complex number (eigenvalue)    |
|                          |                                  | for requested coherent state           |
+--------------------------+----------------------------------+----------------------------------------+
| Coherent density matrix  | ``coherent_dm(N,alpha)``         | same as coherent(N,alpha)              |
| (outer product)          |                                  |                                        |
+--------------------------+----------------------------------+----------------------------------------+
| Thermal density matrix   | ``thermal_dm(N,n)``              | n = particle number expectation value  |
| (for n particles)        |                                  |                                        |
+--------------------------+----------------------------------+----------------------------------------+

and operators:

.. cssclass:: table-striped

+--------------------------+----------------------------+----------------------------------------+
| Operators                | Command (# means optional) | Inputs                                 |
+==========================+============================+========================================+
| Charge operator          | ``charge(N,M=-N)``         | Diagonal operator with entries         |
|                          |                            | from M..0..N.                          |
+--------------------------+----------------------------+----------------------------------------+
| Commutator               | ``commutator(A, B, kind)`` | Kind = 'normal' or 'anti'.             |
+--------------------------+----------------------------+----------------------------------------+
| Diagonals operator       | ``qdiags(N)``              | Quantum object created from arrays of  |
|                          |                            | diagonals at given offsets.            |
+--------------------------+----------------------------+----------------------------------------+
| Displacement operator    | ``displace(N,alpha)``      | N=number of levels in Hilbert space,   |
| (Single-mode)            |                            | alpha = complex displacement amplitude.|
+--------------------------+----------------------------+----------------------------------------+
| Higher spin operators    | ``jmat(j,#s)``             | j = integer or half-integer            |
|                          |                            | representing spin, s = 'x', 'y', 'z',  |
|                          |                            | '+', or '-'                            |
+--------------------------+----------------------------+----------------------------------------+
| Identity                 | ``qeye(N)``                | N = number of levels in Hilbert space. |
+--------------------------+----------------------------+----------------------------------------+
| Lowering (destruction)   | ``destroy(N)``             | same as above                          |
| operator                 |                            |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Momentum operator        | ``momentum(N)``            | same as above                          |
+--------------------------+----------------------------+----------------------------------------+
| Number operator          | ``num(N)``                 | same as above                          |
+--------------------------+----------------------------+----------------------------------------+
| Phase operator           | ``phase(N, phi0)``         | Single-mode Pegg-Barnett phase         |
| (Single-mode)            |                            | operator with ref phase phi0.          |
+--------------------------+----------------------------+----------------------------------------+
| Position operator        | ``position(N)``            | same as above                          |
+--------------------------+----------------------------+----------------------------------------+
| Raising (creation)       | ``create(N)``              | same as above                          |
| operator                 |                            |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Squeezing operator       | ``squeeze(N, sp)``         | N=number of levels in Hilbert space,   |
| (Single-mode)            |                            | sp = squeezing parameter.              |
+--------------------------+----------------------------+----------------------------------------+
| Squeezing operator       | ``squeezing(q1, q2, sp)``  | q1,q2 = Quantum operators (Qobj)       |
| (Generalized)            |                            | sp = squeezing parameter.              |
+--------------------------+----------------------------+----------------------------------------+
| Sigma-X                  | ``sigmax()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Sigma-Y                  | ``sigmay()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Sigma-Z                  | ``sigmaz()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Sigma plus               | ``sigmap()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Sigma minus              | ``sigmam()``               |                                        |
+--------------------------+----------------------------+----------------------------------------+
| Tunneling operator       | ``tunneling(N,m)``         | Tunneling operator with elements of the|
|                          |                            | form :math:`|N><N+m| + |N+m><N|`.      |
+--------------------------+----------------------------+----------------------------------------+

As an example, we give the output for a few of these functions:

.. ipython::

   In [1]: basis(5,3)

   In [2]: coherent(5,0.5-0.5j)

   In [3]: destroy(4)

   In [4]: sigmaz()

   In [5]: jmat(5/2.0,'+')

.. _basics-qobj-props:

Qobj attributes
---------------

We have seen that a quantum object has several internal attributes, such as data, dims, and shape.  These can be accessed in the following way:

.. ipython::

   In [1]: q = destroy(4)

   In [2]: q.dims

   In [3]: q.shape

In general, the attributes (properties) of a ``Qobj`` object (or any Python class) can be retrieved using the `Q.attribute` notation.  In addition to the attributes shown with the ``print`` function, the ``Qobj`` class also has the following:

.. cssclass:: table-striped

+---------------+---------------+----------------------------------------+
| Property      | Attribute     | Description                            |
+===============+===============+========================================+
| Data          | ``Q.data``    | Matrix representing state or operator  |
+---------------+---------------+----------------------------------------+
| Dimensions    | ``Q.dims``    | List keeping track of shapes for       |
|               |               | individual components of a             |
|               |               | multipartite system (for tensor        |
|               |               | products and partial traces).          |
+---------------+---------------+----------------------------------------+
| Shape         | ``Q.shape``   | Dimensions of underlying data matrix.  |
+---------------+---------------+----------------------------------------+
| is Hermitian? | ``Q.isherm``  | Is the operator Hermitian or not?      |
+---------------+---------------+----------------------------------------+
| Type          | ``Q.type``    | Is object of type 'ket, 'bra',         |
|               |               | 'oper', or 'super'?                    |
+---------------+---------------+----------------------------------------+

.. figure:: quide-basics-qobj-box.png
   :align: center
   :width: 3.5in

   The ``Qobj`` Class viewed as a container for the properties need to characterize a quantum operator or state vector.


For the destruction operator above:

.. ipython::

    In [1]: q.type

    In [2]: q.isherm

    In [3]: q.data


The data attribute returns a message stating that the data is a sparse matrix. All ``Qobj`` instances store their data as a sparse matrix to save memory. To access the underlying dense matrix one needs to use the :func:`qutip.Qobj.full` function as described below.

.. _basics-qobj-math:

Qobj Math
----------

The rules for mathematical operations on ``Qobj`` instances are similar to standard matrix arithmetic:

.. ipython::

    In [1]: q = destroy(4)

    In [2]: x = sigmax()

    In [3]: q + 5

    In [4]: x * x

    In [5]: q ** 3

    In [6]: x / np.sqrt(2)


Of course, like matrices, multiplying two objects of incompatible shape throws an error:

.. ipython::
    :okexcept:

    In [1]: q * x


In addition, the logic operators is equal `==` and is not equal `!=` are also supported.

.. _basics-functions:

Functions operating on Qobj class
==================================

Like attributes, the quantum object class has defined functions (methods) that operate on ``Qobj`` class instances. For a general quantum object ``Q``:

.. cssclass:: table-striped

+-----------------+-------------------------------+----------------------------------------+
| Function        | Command                       | Description                            |
+=================+===============================+========================================+
| Check Hermicity | ``Q.check_herm()``            | Check if quantum object is Hermitian   |
+-----------------+-------------------------------+----------------------------------------+
| Conjugate       | ``Q.conj()``                  | Conjugate of quantum object.           |
+-----------------+-------------------------------+----------------------------------------+
| Cosine          | ``Q.cosm()``                  | Cosine of quantum object.              |
+-----------------+-------------------------------+----------------------------------------+
| Dagger (adjoint)| ``Q.dag()``                   | Returns adjoint (dagger) of object.    |
+-----------------+-------------------------------+----------------------------------------+
| Diagonal        | ``Q.diag()``                  | Returns the diagonal elements.         |
+-----------------+-------------------------------+----------------------------------------+
| Diamond Norm    | ``Q.dnorm()``                 | Returns the diamond norm.              |
+-----------------+-------------------------------+----------------------------------------+
| Eigenenergies   | ``Q.eigenenergies()``         | Eigenenergies (values) of operator.    |
+-----------------+-------------------------------+----------------------------------------+
| Eigenstates     | ``Q.eigenstates()``           | Returns eigenvalues and eigenvectors.  |
+-----------------+-------------------------------+----------------------------------------+
| Eliminate States| ``Q.eliminate_states(inds)``  | Returns quantum object with states in  |
|                 |                               | list inds removed.                     |
+-----------------+-------------------------------+----------------------------------------+
| Exponential     | ``Q.expm()``                  | Matrix exponential of operator.        |
+-----------------+-------------------------------+----------------------------------------+
| Extract States  | ``Q.extract_states(inds)``    | Qobj with states listed in inds only.  |
+-----------------+-------------------------------+----------------------------------------+
| Full            | ``Q.full()``                  | Returns full (not sparse) array of     |
|                 |                               | Q's data.                              |
+-----------------+-------------------------------+----------------------------------------+
| Groundstate     | ``Q.groundstate()``           | Eigenval & eigket of Qobj groundstate. |
+-----------------+-------------------------------+----------------------------------------+
| Matrix Element  | ``Q.matrix_element(bra,ket)`` | Matrix element <bra|Q|ket>             |
+-----------------+-------------------------------+----------------------------------------+
| Norm            | ``Q.norm()``                  | Returns L2 norm for states,            |
|                 |                               | trace norm for operators.              |
+-----------------+-------------------------------+----------------------------------------+
| Overlap         | ``Q.overlap(state)``          | Overlap between current Qobj and a     |
|                 |                               | given state.                           |
+-----------------+-------------------------------+----------------------------------------+
| Partial Trace   | ``Q.ptrace(sel)``             | Partial trace returning components     |
|                 |                               | selected using 'sel' parameter.        |
+-----------------+-------------------------------+----------------------------------------+
| Permute         | ``Q.permute(order)``          | Permutes the tensor structure of a     |
|                 |                               | composite object in the given order.   |
+-----------------+-------------------------------+----------------------------------------+
| Projector       | ``Q.proj()``                  | Form projector operator from given     |
|                 |                               | ket or bra vector.                     |
+-----------------+-------------------------------+----------------------------------------+
| Sine            | ``Q.sinm()``                  | Sine of quantum operator.              |
+-----------------+-------------------------------+----------------------------------------+
| Sqrt            | ``Q.sqrtm()``                 | Matrix sqrt of operator.               |
+-----------------+-------------------------------+----------------------------------------+
| Tidyup          | ``Q.tidyup()``                | Removes small elements from Qobj.      |
+-----------------+-------------------------------+----------------------------------------+
| Trace           | ``Q.tr()``                    | Returns trace of quantum object.       |
+-----------------+-------------------------------+----------------------------------------+
| Transform       | ``Q.transform(inpt)``         | A basis transformation defined by      |
|                 |                               | matrix or list of kets 'inpt' .        |
+-----------------+-------------------------------+----------------------------------------+
| Transpose       | ``Q.trans()``                 | Transpose of quantum object.           |
+-----------------+-------------------------------+----------------------------------------+
| Truncate Neg    | ``Q.trunc_neg()``             | Truncates negative eigenvalues         |
+-----------------+-------------------------------+----------------------------------------+
| Unit            | ``Q.unit()``                  | Returns normalized (unit)              |
|                 |                               | vector Q/Q.norm().                     |
+-----------------+-------------------------------+----------------------------------------+

.. ipython::

    In [1]: basis(5, 3)

    In [2]: basis(5, 3).dag()

    In [3]: coherent_dm(5, 1)

    In [4]: coherent_dm(5, 1).diag()

    In [5]: coherent_dm(5, 1).full()

    In [6]: coherent_dm(5, 1).norm()

    In [7]: coherent_dm(5, 1).sqrtm()

    In [8]: coherent_dm(5, 1).tr()

    In [9]: (basis(4, 2) + basis(4, 1)).unit()
