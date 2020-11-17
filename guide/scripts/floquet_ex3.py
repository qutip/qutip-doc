from qutip import *
from scipy import *
import numpy as np
import matplotlib.pyplot as plt

delta = 0.0  * 2*pi; eps0  = 1.0 * 2*pi
A     = 0.15 * 2*pi; omega = 1.0 * 2*pi
T      = (2*pi)/omega
tlist  = linspace(0.0, 20 * T, 501)
psi0   = basis(2,0)

H0 = - delta/2.0 * sigmax() - eps0/2.0 * sigmaz()
H1 = A/2.0 * sigmax()
args = {'w': omega}
H = [H0, [H1, lambda t,args: sin(args['w'] * t)]]

# Noise power spectrum
gamma1 = 0.05
def noise_spectrum(omega):
    if omega>0:
        return 0.5 * gamma1 * omega/(2*pi)
    else:
        return 0

# Collapse operators
   # for fmmesolve
X = sigmax() 

   # for mesolve
c_ops1 = []
gamma = np.zeros([2,2],dtype=complex)
ep, vp = H0.eigenstates()

for i in range(2):
    for j in range(2):
        if i != j:
            gamma[i][j] = 2*np.pi*X.matrix_element(vp[j], vp[i]) \
                          *X.matrix_element(vp[i], vp[j]) \
                          *noise_spectrum(ep[j]-ep[i])
for i in range(2):
    for j in range(2):
        c_ops1.append(sqrt(gamma[i][j])*(vp[i]*vp[j].dag()))

# Find the floquet modes for the time-dependent hamiltonian
f_modes_0, f_energies = floquet_modes(H, T, args)

# Precalculate mode table
f_modes_table_t = floquet_modes_table(f_modes_0, f_energies,
                                      linspace(0, T, 500 + 1), H, T, args)

# Solve the floquet-markov master equation
output = fmmesolve(H, psi0, tlist, [sigmax()], [], [noise_spectrum], T, args)

# Calculate expectation values in the computational basis
p_ex = zeros(shape(tlist), dtype=complex)
for idx, t in enumerate(tlist):
    f_modes_t  = floquet_modes_t_lookup(f_modes_table_t, t, T)
    f_states_t = [np.exp(-1j*t*f_energies[0])*f_modes_t[0], \
                  np.exp(-1j*t*f_energies[1])*f_modes_t[1]]
    p_ex[idx]  = expect(num(2), output.states[idx].transform(f_states_t, True))

# For reference: calculate the same thing with mesolve
output = mesolve(H, psi0, tlist, c_ops1, [num(2)], args)
p_ex_ref = output.expect[0]

# plot the results
from pylab import *
plot(tlist, real(p_ex), 'r--', tlist, 1-real(p_ex), 'b--')
plot(tlist, real(p_ex_ref), 'r', tlist, 1-real(p_ex_ref), 'b')
xlabel('Time')
ylabel('Occupation probability')
legend(("Floquet $P_1$", "Floquet $P_0$", "Lindblad $P_1$", "Lindblad $P_0$"))
show()