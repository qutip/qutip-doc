"""
Basic use of Processor
=============================
 
This example contains the basic functions of
:class:`qutip.qip.device.Processor.`
"""
import copy
import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
from qutip.qip.device import Processor, RandomNoise
from qutip.operators import sigmaz, sigmay, sigmax, destroy
from qutip.metrics import fidelity
from qutip.bloch import Bloch
from qutip.states import basis
from qutip.qip import rx, ry, rz, hadamard_transform

processor = Processor(N=1)
processor.add_ctrl(sigmaz(), targets=0)
processor.add_ctrl(sigmay(), targets=0)
[ctrl for ctrl in processor.ctrls]
processor.coeffs = np.array([[ 0.5, 0.,  0.5],
                            [ 0. , 0.5, 0. ]])
processor.tlist = np.array([0., pi/2., 2*pi/2, 3*pi/2])
result = processor.run_state(rho0=basis(2, 0))
result.states[-1].tidyup(1.0e-5)


tlist = np.linspace(0., 2*np.pi, 20)
processor = Processor(N=1, spline_kind="step_func")
processor.add_ctrl(sigmaz())
processor.tlist = tlist
processor.coeffs = np.array([[np.sin(t) for t in tlist]])
processor.plot_pulses(noisy=False)

tlist = np.linspace(0., 2*np.pi, 20)
processor = Processor(N=1, spline_kind="cubic")
processor.add_ctrl(sigmaz())
processor.tlist = tlist
processor.coeffs = np.array([[np.sin(t) for t in tlist]])
processor.plot_pulses(noisy=False)