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

   In [1]: import matplotlib.pyplot as plt

   In [1]: from qutip.qip.models.spinchain import *

   In [1]: from qutip.qip.models.circuitprocessor import *

   In [1]: from qutip.qip.models.optpulseprocessor import *

Introduction
============
TBA

Qauntum Circuit
=============
TBA

Quantum Circuit Simulator
=========================

In the ``qip`` module, :class:`qutip.CircuitProcessor` is used to simulate quantum circuits at the level of driving Hamiltonian. One could consider it as a simulation of the physical hardware, on which the quantum circuit is to be implemented. As for a physical platform, a processor is determined by a set of Hamiltonians, i.e. the control pulse driving the evolution. Given the control pulse matrix and the corresponding time slices, the evolution can be calculated using the open system solver in QuTiP. The control pulse matrix is a 2D number array, each row is the pulse amplitude for one Hamiltonian and has the same length as the time list. Together they determine the evolution of qubits. One example of defining the control pulse matrix and the time list is as follows:

.. ipython:: 

    In [1]: proc = CircuitProcessor(1)

    In [1]: proc.add_ctrl(sigmaz(), expand_type="periodic")

    In [1]: proc.amps = np.array([[1.0, 1.5, 2.0],
       ...:                       [1.7, 1.3, 0.9]])

    In [1]: proc.tlist = np.array([0.1, 0.3, 0.4])

    In [1]: fig, ax = proc.plot_pulses()

So far, the circuit processor is only a wrap for the open system solver. The key step is how to find the pulse strength and the time duration, under which the evolution is the same as the desired quantum circuit. This gives birth to different subclasses:

| CircuitProcessor
| ├── ModelProcessor
| │   ├── DispersivecQED
| │   └── SpinChain
| └── OptPulseProcessor

In general, there are two ways to find the control pulses. The first one is commonly used in experiments, where a circuit is decomposed into elementary gates, of which the control pulse is well known and can be concatenated to realize the whole quantum circuits. Since the implementation differs from hardware to hardware, it is defined under the class :class:`ModelProcessor` and the corresponding subclasses. The other one is using an algorithm to find the control pulse that leads to the desired evolution. The implementation here makes use of the optimal control module in QuTiP(see :ref:`control`) and is called :class:`OptPulseProcessor`.

Despite this difference, the logic behind circuit processors is the same:

* One defines a circuit processor by a list of available Hamiltonians and, as explained later, hardware dependent noise. The Hamiltonians for a specific model is predefined in the corresponding subclass.

* The control pulse matrix and time slices are either specified by the user or by the method ``load_circuit``, which takes a :class:`qutip.QubitCircuit` or a list of :class:`qutip.Qobj` and will find the control pulse for this circuit.

* The processor calculates the evolution, either analytically or with a numerical solver. In the latter case, collapse operator can be added to simulate decoherence.

If a numerical method is chosen, ``run_state`` in the circuitprocessor returns a object :class:`qutip.solver.Result`. In the case of analytical calculation, a list of the propagator is returned.

SpinChain
---------
TBA

DispersivecQED
--------------
TBA

OptPulseProcessor
-----------------
The :class:`qutip.OptPulseProcessor' uses the method in ``qutip.optimize_pulse_unitary`` in the optimal control module to find the control pulse matrix. All parameters for ``qutip.optimize_pulse_unitary`` can also be given as keyword arguments to ``load_circuit``.

.. ipython:: 

    In [1]: H_d = sigmaz()
       ...: H_c = sigmax()

    In [1]: optproc = OptPulseProcessor(N=1, drift=H_d, ctrls=[H_c])

    In [1]: qc = QubitCircuit(1)
       ...: qc.add_gate("SNOT", targets=[0])

    In [1]: optproc.load_circuit(qc, n_ts=10, evo_time=10)
       ...: real_rho1 = optproc.run_state(rho0=basis(2,0)).states[-1]

    In [1]: ideal_rho1 = (basis(2, 0)+basis(2, 1)).unit()
       ...: fidelity(real_rho1, ideal_rho1)

Noise Simulation
================

The circuit processor in QuTiP, as a circuit simulator, is different from the common simulator of quantum information processing, as it simulates the dynamics of the qubits under the driving Hamiltonian. In many other simulators, the circuit is represented by a sequence of matrices and the noise is simulated with the probabilistic occurrence of collapse. In this sense, circuit processor offers a more physical way of simulating the dynamics in a quantum circuit.

As the design of circuit processor follows the physical realization, so is the noise simulation. Noise can be added to the processor at different levels:

* The decoherence time T1 and T2 can be defined for the processor or for each qubit. When calculating the evolution, the corresponding collapse operators will be added automatically in the solver.

* The noise of the control pulse (e.g. detuned frequency) can be simulated by changing the parameters in the model, e.g. laser frequency in cavity QED.

* It is possible to add further decoherence into the evolution when calling ``run_state`` by adding collapse operators as keyword arguments to the solver.

(There is no real T1 T2 noise in the following code for now since the processor cannot deal with tlist that is not equidistant.)

.. ipython:: 

    In [1]: spinproc = CircularSpinChain(2, T1=0.1, T2=0.5)

    In [1]: qc = QubitCircuit(2)
       ...: qc.add_gate("CNOT", targets=[0], controls=[1])

    In [1]: spinproc.load_circuit(qc)

    In [1]: U_list = spinproc.run_state(states=tensor([basis(2,0),basis(2,1)]))

    In [1]: real_rho = gate_sequence_product(U_list)

    In [1]: ideal_rho = tensor([basis(2,1),basis(2,1)])

    In [1]: fidelity(real_rho, ideal_rho)
