.. QuTiP 
   Copyright (C) 2011-2012, Paul D. Nation & Robert J. Johansson

.. _qip:

*********************************************
Quantum Information Processing
*********************************************

Introduction
============

The Quantum Information Processing (QIP) module aims at providing basic tools for quantum computing simulation both for simple quantum algorithm design and for experiment realization. It offers two different approaches, one with :class:`qutip.qip.QubitCircuit` calculating unitary evolution under quantum gates by matrix production, another called :class:`qutip.qip.CircuitProcessor` using open system solver in QuTiP to simulate NISQ device.

Quantum Circuit
===============

The most common model for quantum computing is the quantum circuit model. In QuTiP, we use :class:`qutip.qip.QubitCircuit` to simulate the unitary evolution. Each quantum gate is saved as a class object :class:`qutip.qip.Gate` with information such as gate name, targets qubits and arguments. To get the matrix representation of each gate, we can call the class method :meth:`qutip.qip.QubitCircuit.propagators()`. Carrying out the matrices production, one gets the result of the evolution or a unitary matrix representation.

.. code-block:: python
   
   >>> from qutip.qip import QubitCircuit, Gate
   >>> qc = QubitCircuit(N=2)
   >>> swap_gate = Gate(name="SWAP", targets=[0, 1])
   >>> qc.add_gate(swap_gate)
   >>> qc.add_gate("CNOT", controls=0, targets=1)
   >>> qc.add_gate(swap_gate)
   >>> print(qc.gates)
   [Gate(SWAP, targets=[0, 1], controls=None), Gate(CNOT, targets=[1],
   controls=[0]), Gate(SWAP, targets=[0, 1], controls=None)]
   >>> U_list = qc.propagators()
   >>> print(gate_sequence_product(U_list))
   Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = True
   Qobj data =
   [[1. 0. 0. 0.]
   [0. 0. 0. 1.]
   [0. 0. 1. 0.]
   [0. 1. 0. 0.]]

The following gates are defined in the class :class:`qutip.qip.Gate`:

====================  ========================================
Gate name                           Description
====================  ========================================
"RX"                  Rotation around x axis
"RY"                  Rotation around y axis
"RZ"                  Rotation around z axis
"SQRTNOT"             Square root of NOT gate
"SNOT"                Hardmard gate
"PHASEGATE"           Add a phase one the state 1
"CRX"                 Controlled rotation around x axis
"CRY"                 Controlled rotation around y axis
"CRZ"                 Controlled rotation around z axis
"CPHASE"              Controlled phase gate
"CNOT"                Controlled NOT gate 
"CSIGN"               Same as CPHASE
"BERKELEY"            Berkeley gate
"SWAPalpha"           SWAPalpha gate
"SWAP"                Swap the states of two qubits
"ISWAP"               Swap gate with additional phase for 01 and 10 states
"SQRTSWAP"            Square root of the SWAP gate
"SQRTISWAP"           Square root of the ISWAP gate
"FREDKIN"             Fredkin gate
"TOFFOLI"             Toffoli gate
"GLOBALPHASE"         Global phase
====================  ========================================

:class:`qutip.qip.QubitCircuit` also has a primitive :meth:`qutip.qip.QubitCircuit.resolve_gates()` method that decomposes the common gates into elementary gate sets such as CNOT or SWAP with single-qubit gates. However, this method is not fully optimized. There is also a method to draw the circuit with LaTeX code. 

In addtion to these pre-defined gates, QuTiP also allow user to define their own gate. The following example shows how to define a customized gate.

.. note::

   Available from QuTiP 4.4

.. code-block:: 

   >>> from qutip.qip import QubitCircuit, Gate, rx
   >>> from qutip import Qobj
   >>> import numpy as np
   >>> def user_gate1(arg_value):
   ...     # controlled rotation X
   ...     mat = np.zeros((4, 4), dtype=np.complex)
   ...     mat[0, 0] = mat[1, 1] = 1.
   ...     mat[2:4, 2:4] = rx(arg_value)
   ...     return Qobj(mat, dims=[[2, 2], [2, 2]])
   ...
   >>> def user_gate2():
   ...     # S gate
   ...     mat = np.array([[1.,   0],
   ...                     [0., 1.j]])
   ...     return Qobj(mat, dims=[[2], [2]])
   ...
   >>>
   >>> qc = QubitCircuit(2)
   >>> qc.user_gates = {"CTRLRX": user_gate1,
   ...                  "S"     : user_gate2}
   ...
   >>> # qubit 0 controlls qubit 1
   ... qc.add_gate("CTRLRX", targets=[0,1], arg_value=np.pi/2)
   >>> # qubit 1 controlls qutbi 0
   ... qc.add_gate("CTRLRX", targets=[1,0], arg_value=np.pi/2)
   >>> # a gate can also be added using the Gate class
   ... g_T = Gate("S", targets=[1])
   >>> qc.add_gate("S", targets=[1])
   >>> props = qc.propagators()
   >>> props[0]
   Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = False
   Qobj data =
   [[1.+0.j 0.+0.j 0.+0.j 0.+0.j]
   [0.+0.j 0.+1.j 0.+0.j 0.+0.j]
   [0.+0.j 0.+0.j 1.+0.j 0.+0.j]
   [0.+0.j 0.+0.j 0.+0.j 0.+1.j]]
   >>> props[1]
   Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = False
   Qobj data =
   [[1.+0.j 0.+0.j 0.+0.j 0.+0.j]
   [0.+0.j 0.+1.j 0.+0.j 0.+0.j]
   [0.+0.j 0.+0.j 1.+0.j 0.+0.j]
   [0.+0.j 0.+0.j 0.+0.j 0.+1.j]]
   >>> props[2]
   Quantum object: dims = [[2, 2], [2, 2]], shape = (4, 4), type = oper, isherm = False
   Qobj data =
   [[1.        +0.j         0.        +0.j         0.        +0.j
   0.        +0.j        ]
   [0.        +0.j         1.        +0.j         0.        +0.j
   0.        +0.j        ]
   [0.        +0.j         0.        +0.j         0.70710678+0.j
   0.        -0.70710678j]
   [0.        +0.j         0.        +0.j         0.        -0.70710678j
   0.70710678+0.j        ]]

CircuitProcessor for QIP simulation
===================================

.. note::

   Available from QuTiP 4.5

In addition to direct matrix production, QuTiP also has anther approach to QIP simulation. Based on the open system solver, :class:`qutip.qip.CircuitProcessor` in the :mod:`qutip.qip` module simulates quantum circuits at the level of driving Hamiltonians. One can consider the circuit processor as a simulator of a quantum device, on which the quantum circuit is to be implemented. Like a real quantum device, the processor is determined by a list of Hamiltonians, i.e. the control pulse driving the evolution. Given the intensity of the control pulses and the corresponding time slices for each pulse, the evolution can be calculated using the solver. The pulse intensity and time for each pulse are saved in the attributes :attr:`qutip.qip.CircuitProcessor.coeffs`, a 2-d NumPy array, and :attr:`qutip.qip.CircuitProcessor.tlist`, a 1-d NumPy array. We can either use the coefficients as a step function or with cubic spline. For step function, tlist specifies the start and the end of each pulse and thus is one element longer the coeffs. One example of defining the control pulse coefficients and the time array is as follows:

.. code-block:: python

   >>> from qutip.qip import CircuitProcessor
   >>> proc = CircuitProcessor(2)
   >>> proc.add_ctrl(sigmaz(), cyclic_permutation=True)  # for all qubits
   >>> proc.coeffs = np.array([[1.0, 1.5, 2.0], [1.8, 1.3, 0.8]])
   >>> proc.tlist = np.array([0.1, 0.2, 0.4, 0.5])

This is the framework and most essential part of the simulator's API, for now, it looks like just a wrap for the open system solver. However, based on this, we can implement different physical realization. They differ mainly in how to find the control pulse for a quantum circuit, which gives birth to different sub-classes:

| CircuitProcessor
| ├── ModelProcessor
| │   ├── DispersivecQED
| │   └── SpinChain
| └── OptPulseProcessor

In general, there are two ways to find the control pulses. The first one is more experiment oriented and based on physical models. We define a universal set of
gates in the processor as well as how to implement them on the physical hardware. This is usually the case where control pulses realizing those gates are well known and can be concatenated to realize the whole quantum circuits. Since they are based on physical models, they are called :class:`qutip.qip.ModelProcessor`. Two realizations have already been implemented: the spin chain and the CQED model for quantum computing. In those models, the driving Hamiltonians are predefined. The other approach, based on the optimal control module in QuTiP (see :ref:`control`), is called :class:`qutip.qip.OptPulseProcessor`. In this framework, one can define the available Hamiltonians in their system. The processor then uses algorithms to find the optimal control pulse that leads to the desired unitary evolution.

Despite this difference, the logic behind circuit processors is all the same:

* One defines a circuit processor by a list of available Hamiltonians and, as explained later, hardware-dependent noise. In model bases processor, the Hamiltonians are predefined and one only need to give the device parameters like frequency and interaction strength. 

* The control pulse coefficients and time slices are either specified by the user or calculated by the method :meth:`qutip.qip.CircuitProcessor.load_circuit()`, which takes a :class:`qutip.qip.QubitCircuit` and find the control pulse for this evolution.

* The processor calculates the evolution, either analytically (without noise) or with a numerical solver. In the latter case, collapse operator can be added to simulate decoherence.

If a numerical method is chosen, :meth:`qutip.qip.CircuitProcessor.run_state` returns a object :class:`qutip.solver.Result`. In the case of analytical calculation, a list of the propagators is returned.

SpinChain
---------

:class:`qutip.qip.LinearSpinChain` and :class:`qutip.qip.CircularSpinChain` are quantum computing models base on the spin chain realization. The control Hamiltonians are :math:`\sigma_x`, :math:`\sigma_z` and :math:`\sigma_x \sigma_x + \sigma_y \sigma_y`. This processor will first decompose the gate into the universal gate set with ISWAP and SQRTISWAP as two-qubits gates, resolve them into quantum gates of adjacent qubits and then calculate the pulse coefficients.

DispersivecQED
--------------

Same as above, :class:`qutip.qip.DispersivecQED` is a representation base on Cavity Quantum Electrodynamics. The workflow is similar to the one for the spin chain, except that the component systems are a multi-level cavity and a qubits system. The control Hamiltonians are the single qubits rotation together with the qubits-cavity interaction :math:`a^{\dagger} \sigma^{-} + a \sigma^{+}`. The device parameters including the cavity frequency, qubits frequency, detuning and interaction strength etc.

OptPulseProcessor
-----------------
The :class:`qutip.qip.OptPulseProcessor` uses the function in :func:`qutip.control.pulseoptim.optimize_pulse_unitary` in the optimal control module to find the control pulse matrix. The Hamiltonian including a drift part and a control part and only the control part will be optimized. The unitary evolution then goes as

.. math::

   U(\Delta t)=\exp(\rm{i} \cdot \Delta t [H_d  + \sum_j u_j H_j] )

All parameters for :func:`qutip.control.pulseoptim.optimize_pulse_unitary` can also be given as keyword arguments to the class method :meth:`qutip.qip.OptPulseProcessor.load_circuit` for advanced requirements.

.. code-block:: python

   >>> from qutip.qip.models.optpulseprocessor import OptPulseProcessor
   >>> H_d = sigmaz()
   >>> H_c = sigmax()
   >>> optproc = OptPulseProcessor(N=1, drift=H_d, ctrls=[H_c])
   >>> qc = QubitCircuit(1)
   >>> qc.add_gate("SNOT", targets=[0])
   >>> optproc.load_circuit(qc, n_ts=10, evo_time=10)
   (array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10.]),
   array([[ 0.30731151, -0.65352563, -0.48110828,  0.28911593, -0.29193726,
         -0.06580244, -0.80886723, -0.31633408, -0.66637145,  0.19453724]]))
   >>> real_rho1 = optproc.run_state(rho0=basis(2,0)).states[-1]
   >>> ideal_rho1 = (basis(2, 0) + basis(2, 1)).unit()
   >>> fidelity(real_rho1, ideal_rho1)

Noise Simulation
================

In the common way of QIP simulation, where evolution is carried out by gate matrix production, the noise is usually simulated with bit flipping and sign flipping errors. The typical approaches are either applying bit/sign flipping gate probabilistically or applying Kraus operators representing different noisy channels (e.g. amplitude damping, dephasing) after each unitary gate evolution. Those two ways are equivalent and the parameters in the Kraus operators are exactly the probability of a flipping error happens during the gate operation time.

Since the circuit processor simulates the state evolution at the level of driving Hamiltonian, there is no way to apply an error operator. Instead, the error is directly added to the driving Hamiltonian list or the collapse operators. Mathematically, this is no different from adding flipping error probabilistically (It is actually how :func:`qutip.mcsolve` works). The collapse operator for amplitude damping and dephasing are exactly the destroying operator and sign-flipping operator. One just needs to choose the correct coefficients for them to simulate, e,g, the relaxation time T1, T2. This simulator is closer to the physical implementation and also more general because it is based on the open system evolution instead of abstract operators.

Compared to the approach of Kraus operators, this way of simulating noise is more computationally expensive. If you only want to simulate the decoherence of single-qubit relaxation, there is no need to go through this approach. However, it is closer to the real experimental and, therefore, more convenient in some cases, such as when the noise interested is not limited to single-qubit relaxation. For instance, a pulse on one qubit might affect the neighbouring qubits, the evolution is still unitary but the gate fidelity will decrease. It is not always easy or even possible to define a noisy gate matrix. In this simulator, it can be easily down by defining a :class:`qutip.qip.ControlAmpNoise`. Here we show two examples:

The first example is a processor with one qubit under rotation around the z-axis with relaxation time T2=5. We can measure the population of the :math:`\left| + \right\rangle` state and observe the Ramsey signal:

.. image:: /gallery/auto_examples/qip/images/sphx_glr_plot_qip_relaxation_001.png

The second example demonstrates a biased Gaussian noise on the pulse amplitude. For visualization purpose, we plot the noisy pulse intensity instead of the state fidelity. The three pulses can, for example, be a zyz-decomposition of an arbitrary single-qubit gate:

.. image:: /gallery/auto_examples/qip/images/sphx_glr_plot_qip_amplitude_noise_001.png

.. image:: /gallery/auto_examples/qip/images/sphx_glr_plot_qip_amplitude_noise_001.png

As the design of circuit processor follows the physical realization, so is the noise simulation. Noise can be added to the processor at different levels:

* The decoherence time T1 and T2 can be defined for the processor or each qubit. When calculating the evolution, the corresponding collapse operators will be added automatically in the solver.

* The noise of the physical parameters (e.g. detuned frequency) can be simulated by changing the parameters in the model, e.g. laser frequency in cavity QED. (This can only be time-independent since QuTiP open system solver only allowing varying coefficients, not varying Hamiltonian operators.)

* The noise of the pulse intensity can be simulated by modifying the coefficients of the Hamiltonian operators or even adding new Hamiltonians.

The simplest relaxation noise can be defined directly in the circuit processor. There are a few predefined noise objects :class:`qutip.qip.CircuitNoise` that can be added to the simulation with the method :meth:`qutip.qip.CircuitProcessor.add_noise`.

Workflow of the CircuitProcessor
================================
If you are interested in the workflow inside the simulator, you can have a look at the two figures bellow.

The figure above shows how the noise is processed in the circuit processor. The noise is defined separately in a class object. When called, it takes parameters and unitary noiseless :class:`qutip.QobjEvo` from the processor, generated the noisy version and send the noise :class:`qutip.QobjEvo` together with the collapse operators to the processor.

.. image:: /figures/qip/CircuitProcessor-workflow.png

When calculating the evolution, the processor first creates its own :class:`qutip.QobjEvo` of the noiseless desired evolution. 
It will then find all the noise objects saved the attribute :attr:`qutip.qip.CircuitProcessor.noise` and call the corresponding methods to get the :class:`qutip.QobjEvo` and a list of collapse operatrors representing these noisy dynamics.(For collapse operators, we don't want to add all the constant collapse into one time-independent operator, so we use a list). 
The processor then combines its own :class:`qutip.QobjEvo` with those from the noise object and give them to the solver.

.. image:: /figures/qip/CircuitProcessor-noise.png