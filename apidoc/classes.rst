.. QuTiP
   Copyright (C) 2011-2012, Paul D. Nation & Robert J. Johansson

.. _classes:

***************
Classes
***************

.. _classes-qobj:

Qobj
--------------

.. autoclass:: qutip.Qobj
    :members:

.. _classes-qobjevo:

QobjEvo
--------------

.. autoclass:: qutip.QobjEvo
    :members:

.. _classes-eseries:

eseries
-----------------

.. autoclass:: qutip.eseries
    :members:

.. _classes-bloch:

Bloch sphere
---------------

.. autoclass:: qutip.bloch.Bloch
    :members:


Cubic Spline
---------------

.. autoclass:: qutip.interpolate.Cubic_Spline
    :members:


.. _classes-non_markov:

Non-Markovian Solvers
---------------------

.. autoclass:: qutip.nonmarkov.heom.HEOMSolver
    :members:

.. autoclass:: qutip.nonmarkov.heom.HSolverDL
    :members:

.. autoclass:: qutip.nonmarkov.memorycascade.MemoryCascade
    :members:

.. autoclass:: qutip.nonmarkov.transfertensor.TTMSolverOptions
    :members:

.. _classes-odeoptions:

Solver Options and Results
---------------------------

.. autoclass:: qutip.solver.Options
    :members:

.. autoclass:: qutip.solver.Result
    :members:

.. autoclass:: qutip.solver.Stats
    :members:

.. autoclass:: qutip.stochastic.StochasticSolverOptions
    :members:

.. _classes-piqs:

Permutational Invariance
---------------

.. autoclass:: qutip.piqs.Dicke
    :members:

.. autoclass:: qutip.piqs.Pim
    :members:

.. _classes-distributions:

Distribution functions
----------------------------

.. autoclass:: qutip.distributions.Distribution
    :members:

.. autoclass:: qutip.distributions.WignerDistribution
    :members:

.. autoclass:: qutip.distributions.QDistribution
    :members:

.. autoclass:: qutip.distributions.TwoModeQuadratureCorrelation
    :members:

.. autoclass:: qutip.distributions.HarmonicOscillatorWaveFunction
    :members:

.. autoclass:: qutip.distributions.HarmonicOscillatorProbabilityFunction
    :members:

.. _classes-qip:

Quantum information processing
------------------------------

.. autoclass:: qutip.qip.Gate
    :members:

.. autoclass:: qutip.qip.QubitCircuit
    :members:

.. autoclass:: qutip.qip.CircuitProcessor
    :members:

.. autoclass:: qutip.qip.OptPulseProcessor
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.ModelProcessor
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.SpinChain
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.LinearSpinChain
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.CircularSpinChain
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.DispersivecQED
    :members:
    :inherited-members:
    
.. autoclass:: qutip.qip.CircuitNoise
    :members:

.. autoclass:: qutip.qip..DecoherenceNoise
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.RelaxationNoise
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.ControlAmpNoise
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.WhiteNoise
    :members:
    :inherited-members:

.. autoclass:: qutip.qip.UserNoise
    :members:
    :inherited-members:

.. _classes-control:

Optimal control
---------------

.. autoclass:: qutip.control.optimizer.Optimizer
    :members:

.. autoclass:: qutip.control.optimizer.OptimizerBFGS
    :members:

.. autoclass:: qutip.control.optimizer.OptimizerLBFGSB
    :members:

.. autoclass:: qutip.control.optimizer.OptimizerCrab
    :members:

.. autoclass:: qutip.control.optimizer.OptimizerCrabFmin
    :members:

.. autoclass:: qutip.control.optimizer.OptimIterSummary
    :members:

.. autoclass:: qutip.control.termcond.TerminationConditions
    :members:

.. autoclass:: qutip.control.optimresult.OptimResult
    :members:

.. autoclass:: qutip.control.dynamics.Dynamics
    :members:

.. autoclass:: qutip.control.dynamics.DynamicsGenMat
    :members:

.. autoclass:: qutip.control.dynamics.DynamicsUnitary
    :members:

.. autoclass:: qutip.control.dynamics.DynamicsSymplectic
    :members:

.. autoclass:: qutip.control.propcomp.PropagatorComputer
    :members:

.. autoclass:: qutip.control.propcomp.PropCompApproxGrad
    :members:

.. autoclass:: qutip.control.propcomp.PropCompDiag
    :members:

.. autoclass:: qutip.control.propcomp.PropCompFrechet
    :members:

.. autoclass:: qutip.control.fidcomp.FidelityComputer
    :members:

.. autoclass:: qutip.control.fidcomp.FidCompUnitary
    :members:

.. autoclass:: qutip.control.fidcomp.FidCompTraceDiff
    :members:

.. autoclass:: qutip.control.fidcomp.FidCompTraceDiffApprox
    :members:

.. autoclass:: qutip.control.tslotcomp.TimeslotComputer
    :members:

.. autoclass:: qutip.control.tslotcomp.TSlotCompUpdateAll
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGen
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenRandom
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenZero
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenLinear
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenPeriodic
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenSine
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenSquare
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenSaw
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenTriangle
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenGaussian
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenGaussianEdge
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenCrab
    :members:

.. autoclass:: qutip.control.pulsegen.PulseGenCrabFourier
    :members:

.. autoclass:: qutip.control.stats.Stats
    :members:

.. autoclass:: qutip.control.dump.Dump
    :members:

.. autoclass:: qutip.control.dump.OptimDump
    :members:

.. autoclass:: qutip.control.dump.DynamicsDump
    :members:

.. autoclass:: qutip.control.dump.DumpItem
    :members:

.. autoclass:: qutip.control.dump.EvoCompDumpItem
    :members:

.. autoclass:: qutip.control.dump.DumpSummaryItem
    :members:
