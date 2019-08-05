"""
Basic use of CircuitProcessor
=============================
 
This example contains the basic function of CircuitProcessor.
"""
import copy
import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
from qutip.qip import CircuitProcessor
from qutip.qip import WhiteNoise
from qutip import sigmaz, sigmay, sigmax, destroy, basis, fidelity, Bloch
from qutip.qip import rx, ry, rz, hadamard_transform

processor = CircuitProcessor(N=1)
processor.add_ctrl(sigmaz(), targets=0)
processor.add_ctrl(sigmay(), targets=0)
[ctrl for ctrl in processor.ctrls]
processor.coeffs = np.array([[ 0.5, 0.,  0.5],
                            [ 0. , 0.5, 0. ]])
processor.tlist = np.array([0., pi/2., 2*pi/2, 3*pi/2])
result = processor.run_state(rho0=basis(2, 0))
result.states[-1].tidyup(1.0e-5)


tlist = np.linspace(0., 2*np.pi, 20)
processor = CircuitProcessor(N=1, spline_kind="step_func")
processor.add_ctrl(sigmaz())
processor.tlist = tlist
processor.coeffs = np.array([[np.sin(t) for t in tlist]])
processor.plot_pulses(noisy=False)

tlist = np.linspace(0., 2*np.pi, 20)
processor = CircuitProcessor(N=1, spline_kind="cubic")
processor.add_ctrl(sigmaz())
processor.tlist = tlist
processor.coeffs = np.array([[np.sin(t) for t in tlist]])
processor.plot_pulses(noisy=False)