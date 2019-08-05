"""
Control Amplitude noise
=======================

This example demonstrate how to add Gaussian noise to the control pulse.
"""
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
processor.plot_pulses(noisy=False, title="Original control amplitude", figsize=(5,3))
processor.plot_pulses(noisy=True, title="Noisy control amplitude", figsize=(5,3))