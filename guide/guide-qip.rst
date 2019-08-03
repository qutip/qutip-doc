.. QuTiP 
   Copyright (C) 2011-2012, Paul D. Nation & Robert J. Johansson

.. _qip:

*********************************************
Quantum Information Processing
*********************************************

.. ipython::
   :suppress:

   In [1]: from qutip import *

   In [1]: import numpy as np

Introduction
============

The Quantum Information Processing (QIP) module aims at providing basic tools for quantum computing simulation both for simple quantum algorithm design and for experiment realization. It offers two different approaches, one with :class:`qutip.qip.QubitCircuit` calculating unitary evolution under quantum gates by matrix production, another called :class:`qutip.qip.CircuitProcessor` using open system solver in QuTiP to simulate NISQ device.

Qauntum Circuit
===============

The most common model for quantum computing is the quantum circuit model. In QuTiP, we use :class:`qutip.qip.QubitCircuit` to simulate the unitary evolution. Each quantum gate is saved as a class object :class:`qutip.qip.Gate` with information such as gate name, targets qubits and arguments. To get the matrix representation of each gate, we can call the class method :meth:`qutip.qip.Circuit.get_propagators`. Carrying out the matrices production, one gets the result of the evolution or a unitary matrix representation.

.. IPython::
   
   In [1]: from qutip.qip import QubitCircuit, Gate

   In [1]: qc = QubitCircuit(N=2)

   In [1]: swap_gate = Gate(name="SWAP", targets=[0, 1])

   In [1]: qc.add_gate(swap_gate)
      ...: qc.add_gate("CNOT", controls=0, targets=1)
      ...: qc.add_gate(swap_gate)

   In [1]: qc.gates

   In [1]: U_list = qc.propagators()

   In [1]: gate_sequence_product(U_list)

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

:class:`qutip.qip.QubitCircuit` also has a primitive :meth:`qutip.qip.QubitCircuit.resolve_gates` method that decomposes the common gates into elementary gate sets such as CNOT or SWAP with single-qubit gates. However, this method is not fully optimized. There is also a method to draw the circuit with LaTeX code. 

CircuitProcessor for QIP simulation
===================================

In addition to direct matrix production, QuTiP also has anther approach to QIP simulation. Based on the open system solver, :class:`qutip.qip.CircuitProcessor` in the :mod:`qutip.qip` module simulates quantum circuits at the level of driving Hamiltonians. One can consider the circuit processor as a simulator of a quantum device, on which the quantum circuit is to be implemented. Like a real quantum device, the processor is determined by a list of Hamiltonians, i.e. the control pulse driving the evolution. Given the intensity of the control pulses and the corresponding time slices for each pulse, the evolution can be calculated using the solver. The pulse intensity and time for each pulse are saved in the attributes :attr:`qutip.qip.CircuitProcessor.coeffs`, a 2-d NumPy array, and :attr:`qutip.qip.CircuitProcessor.tlist`, a 1-d NumPy array. We can either use the coefficients as a step function or with cubic spline. For step function, tlist specifies the start and the end of each pulse and thus is one element longer the coeffs. One example of defining the control pulse coefficients and the time array is as follows:

.. ipython::

    In [1]: from qutip.qip import CircuitProcessor

    In [1]: proc = CircuitProcessor(2)

    In [1]: proc.add_ctrl(sigmaz(), cyclic_permutation=True)  # for all qubits

    In [1]: proc.coeffs = np.array([[1.0, 1.5, 2.0], [1.8, 1.3, 0.8]])

    In [1]: proc.tlist = np.array([0.1, 0.2, 0.4, 0.5])

.. plot::

   from qutip import *
   import matplotlib.pyplot as plt
   import numpy as np
   from qutip.qip import CircuitProcessor
   proc = CircuitProcessor(2)
   proc.add_ctrl(sigmaz(), cyclic_permutation=True)
   proc.coeffs = np.array([[1.0, 1.5, 2.0], [1.8, 1.3, 0.8]])
   proc.tlist = np.array([0.1, 0.2, 0.4, 0.5])
   proc.plot_pulses(title="Control pulses")
   plt.show()

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

* The control pulse coefficients and time slices are either specified by the user or calculated by the method :meth:`qutip.qip.CircuitProcessor.load_circuit`, which takes a :class:`qutip.qip.QubitCircuit` and find the control pulse for this evolution.

* The processor calculates the evolution, either analytically (without noise) or with a numerical solver. In the latter case, collapse operator can be added to simulate decoherence.

If a numerical method is chosen, :meth:`qutip.qip.CircuitProcessor.run_state` returns a object :class:`qutip.solver.Result`. In the case of analytical calculation, a list of the propagators is returned.

SpinChain
---------

:class:`qutip.qip.LinearSpinChain` and :class:`qutip.qip.CircularSpinChain` are quantum computing models base on the spin chain realization. The control Hamiltonians are :math:`\sigma_x`, :math:`\sigma_z` and :math:`\sigma_x \sigma_x + \sigma_y \sigma_y`. This processor will first decompose the gate into the universal gate set with ISWAP and SQRTISWAP as two-qubits gates, resolve them into quantum gates of adjacent qubits and then calculate the pulse coefficients.

DispersivecQED
--------------

Same as above, :class:`qutip.DispersivecQED` is a representation base on Cavity Quantum Electrodynamics. The workflow is similar to the one for the spin chain, except that the component systems are a multi-level cavity and a qubits system. The control Hamiltonians are the single qubits rotation together with the qubits-cavity interaction :math:`a^{\dagger} \sigma^{-} + a \sigma^{+}`. The device parameters including the cavity frequency, qubits frequency, detuning and interaction strength etc.

OptPulseProcessor
-----------------
The :class:`qutip.OptPulseProcessor` uses the function in :func:`qutip.control.optimize_pulse_unitary` in the optimal control module to find the control pulse matrix. The Hamiltonian including a drift part and a control part and only the control part will be optimized. The unitary evolution then goes as

.. math::

   U(\Delta t)=\exp(\rm{i} \cdot \Delta t [H_d  + \sum_j u_j H_j] )

All parameters for :func:`qutip.control.optimize_pulse_unitary` can also be given as keyword arguments to the class method :meth:`qutip.OptPulseProcessor.load_circuit` for advanced requirements.

.. ipython:: 

    In [1]: from qutip.qip.models.optpulseprocessor import OptPulseProcessor

    In [1]: H_d = sigmaz()
       ...: H_c = sigmax()

    In [1]: optproc = OptPulseProcessor(N=1, drift=H_d, ctrls=[H_c])

    In [1]: qc = QubitCircuit(1)
       ...: qc.add_gate("SNOT", targets=[0])

    In [1]: optproc.load_circuit(qc, n_ts=10, evo_time=10)
       ...: real_rho1 = optproc.run_state(rho0=basis(2,0)).states[-1]

    In [1]: ideal_rho1 = (basis(2, 0) + basis(2, 1)).unit()
       ...: fidelity(real_rho1, ideal_rho1)

Noise Simulation
================

In the common way of QIP simulation, where evolution is carried out by gate matrix production, the noise is usually simulated with bit flipping and sign flipping errors. The typical approaches are either applying bit/sign flipping gate probabilistically or applying Kraus operators representing different noisy channels (e.g. amplitude damping, dephasing) after each unitary gate evolution. Those two ways are equivalent and the parameters in the Kraus operators are exactly the probability of a flipping error happens during the gate operation time.

Since the circuit processor simulates the state evolution at the level of driving Hamiltonian, there is no way to apply an error operator. Instead, the error is directly added to the driving Hamiltonian list or the collapse operators. Mathematically, this is no different from adding flipping error probabilistically (It is actually how :func:`qutip.mcsolve` works). The collapse operator for amplitude damping and dephasing are exactly the destroying operator and sign-flipping operator. One just needs to choose the correct coefficients for them to simulate, e,g, the relaxation time T1, T2. This simulator is closer to the physical implementation and also more general because it is based on the open system evolution instead of abstract operators.

Compared to the approach of Kraus operators, this way of simulating noise is more computationally expensive. If you only want to simulate the decoherence of single-qubit relaxation, there is no need to go through this approach. However, it is closer to the real experimental and, therefore, more convenient in some cases, such as when the noise interested is not limited to single-qubit relaxation. For instance, a pulse on one qubit might affect the neighbouring qubits, the evolution is still unitary but the gate fidelity will decrease. It is not always easy or even possible to define a noisy gate matrix. In this simulator, it can be easily down by defining a :class:`qutip.qip.ControlAmplitudeNoise`. Here we show two examples:

The first example is a processor with one qubit under rotation around the z-axis with relaxation time T2=5. We can measure the population of the :math:`\left| + \right\rangle` state and observe the Ramsey signal:

.. plot::

   import numpy as np
   import matplotlib.pyplot as plt
   from qutip.qip import CircuitProcessor
   from qutip.operators import sigmaz, destroy
   from qutip.qip import snot
   from qutip.states import basis
   a = destroy(2)
   Hadamard = snot()
   plus_state = (basis(2,1) + basis(2,0)).unit()
   tlist = np.arange(0.00, 20.2, 0.2)

   T2 = 5
   processor = CircuitProcessor(1, T2=T2)
   processor.add_ctrl(sigmaz())
   processor.tlist= tlist
   processor.coeffs = np.ones((1,len(processor.tlist)))
   result = processor.run_state(
      rho0=plus_state, e_ops=[a.dag()*a, Hadamard*a.dag()*a*Hadamard])

   fig, ax = plt.subplots()
   # detail about length of tlist needs to be fixed
   ax.plot(tlist[:-1], result.expect[1][:-1], '.', label="simulation")
   ax.plot(tlist[:-1], np.exp(-1./T2*tlist[:-1])*0.5 + 0.5, label="theory")
   ax.set_xlabel("t")
   ax.set_ylabel("Ramsey signal")
   ax.legend()
   ax.set_title("Relaxation T2=5")
   ax.grid()
   fig.tight_layout()
   plt.show()

The second example demonstrates a biased Gaussian noise on the pulse amplitude. For visualization purpose, we plot the noisy pulse intensity instead of the state fidelity. The three pulses can, for example, be a zyz-decomposition of an arbitrary single-qubit gate:

.. plot::

   import numpy as np
   import matplotlib.pyplot as plt
   from qutip.operators import sigmaz, sigmay
   from qutip.qip import CircuitProcessor
   from qutip.qip import WhiteNoise

   # add control Hamiltonians
   processor = CircuitProcessor(N=1)
   processor.add_ctrl(sigmaz(), targets=0)
   processor.add_ctrl(sigmay(), targets=0)

   # define coeffs and tlist
   processor.coeffs = np.array([[ 0.3, 0.,  0.2],
                            [ 0. , 0.5, 0. ]])
   processor.tlist = np.array([0., np.pi/2., 2*np.pi/2, 3*np.pi/2])

   # define noise
   processor.add_noise(WhiteNoise(mean=0.08, std=0.02, dt=0.1))
   processor.plot_pulses(noisy=False, title="Original control amplitude")
   processor.plot_pulses(noisy=True, title="Noisy control amplitude")
   plt.show()

As the design of circuit processor follows the physical realization, so is the noise simulation. Noise can be added to the processor at different levels:

* The decoherence time T1 and T2 can be defined for the processor or each qubit. When calculating the evolution, the corresponding collapse operators will be added automatically in the solver.

* The noise of the physical parameters (e.g. detuned frequency) can be simulated by changing the parameters in the model, e.g. laser frequency in cavity QED. (This can only be time-independent since QuTiP open system solver only allowing varying coefficients, not varying Hamiltonian operators.)

* The noise of the pulse intensity can be simulated by modifying the coefficients of the Hamiltonian operators or even adding new Hamiltonians.

The simplest relaxation noise can be defined directly in the circuit processor. There are a few predefined noise objects :class:`qutip.qip.CircuitNoise` that can be added to the simulation with the method :meth:`qutip.qip.CircuitProcessor.add_noise`.
