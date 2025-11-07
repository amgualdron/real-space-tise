import numpy as np
import matplotlib.pyplot as plt

# Load the wavefunction data
# Skip the header lines (starting with #)
data = np.loadtxt('results/qho_wavefunctions.txt', comments='#')

# Extract columns
x = data[:, 0]      # Position
psi1 = data[:, 1]   # Ground state
psi2 = data[:, 2]   # First excited state
psi3 = data[:, 3]   # Second excited state
psi4 = data[:, 4]   # Third excited state
psi5 = data[:, 5]   # Fourth excited state

# Calculate probability densities
prob1 = psi1**2
prob2 = psi2**2
prob3 = psi3**2
prob4 = psi4**2
prob5 = psi5**2

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Plot 1: Wavefunctions
ax1.plot(x, psi1, label='$\psi_0$ (n=0)', linewidth=2)
ax1.plot(x, psi2, label='$\psi_1$ (n=1)', linewidth=2)
ax1.plot(x, psi3, label='$\psi_2$ (n=2)', linewidth=2)
ax1.plot(x, psi4, label='$\psi_3$ (n=3)', linewidth=2)
ax1.plot(x, psi5, label='$\psi_4$ (n=4)', linewidth=2)
ax1.set_xlabel('Position x (Å)', fontsize=12)
ax1.set_ylabel('Wavefunction $\psi(x)$', fontsize=12)
ax1.set_title('Quantum Harmonic Oscillator Eigenfunctions', fontsize=14)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim([x.min(), x.max()])

# Plot 2: Probability Densities
ax2.plot(x, prob1, label='$|\psi_0|^2$ (n=0)', linewidth=2)
ax2.plot(x, prob2, label='$|\psi_1|^2$ (n=1)', linewidth=2)
ax2.plot(x, prob3, label='$|\psi_2|^2$ (n=2)', linewidth=2)
ax2.plot(x, prob4, label='$|\psi_3|^2$ (n=3)', linewidth=2)
ax2.plot(x, prob5, label='$|\psi_4|^2$ (n=4)', linewidth=2)
ax2.set_xlabel('Position x (Å)', fontsize=12)
ax2.set_ylabel('Probability Density $|\psi(x)|^2$, fontsize=12')
ax2.set_title('Probability Densities', fontsize=14)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim([x.min(), x.max()])

# Optional: zoom into the central region where wavefunctions are significant
# ax1.set_xlim([-10, 10])
# ax2.set_xlim([-10, 10])

plt.tight_layout()
plt.savefig('docs/figures/qho_wavefunctions.png', dpi=300)
plt.show()

print("Plot saved as 'results/qho_wavefunctions.png'")