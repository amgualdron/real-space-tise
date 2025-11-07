import numpy as np
import matplotlib.pyplot as plt

# Load the wavefunction data
data = np.loadtxt('results/qho_wavefunctions.txt', comments='#')

# Extract columns
x = data[:, 0]
psi1 = data[:, 1]
psi2 = data[:, 2]
psi3 = data[:, 3]
psi4 = data[:, 4]
psi5 = data[:, 5]

# Calculate probability densities
prob1 = psi1**2
prob2 = psi2**2
prob3 = psi3**2
prob4 = psi4**2
prob5 = psi5**2

# Constants (matching your Fortran code)
hbar_omega = 1.0  # eV
hbar = 6.582119569e-16  # eV·s
c = 2.997924580e18  # Å/s
m_e_c2 = 0.51099895e6  # eV
m_e = m_e_c2 / c**2  # eV·s²/Ų
omega = hbar_omega / hbar  # rad/s

V = 0.5 * m_e * omega**2 * x**2  # Potential in eV

# Energy levels (from your results)
E = np.array([0.499990, 1.499949, 2.499867, 3.499743, 4.499579])

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Wavefunctions with potential and energy levels
ax1.plot(x, psi1, label='$\psi_0$ (n=0)', linewidth=2)
ax1.plot(x, psi2, label='$\psi_1$ (n=1)', linewidth=2)
ax1.plot(x, psi3, label='$\psi_2$ (n=2)', linewidth=2)
ax1.plot(x, psi4, label='$\psi_3$ (n=3)', linewidth=2)
ax1.plot(x, psi5, label='$\psi_4$ (n=4)', linewidth=2)

# Add potential on secondary y-axis
ax1_right = ax1.twinx()
ax1_right.plot(x, V, 'k--', linewidth=2.5, alpha=0.6, label='V(x)')
ax1_right.set_ylabel('Potential V(x) (eV)', fontsize=12)

# Add horizontal lines for energy levels
colors = ['C0', 'C1', 'C2', 'C3', 'C4']
for i, (energy, color) in enumerate(zip(E, colors)):
    ax1_right.axhline(y=energy, color=color, linestyle=':', linewidth=1.5, alpha=0.7)

# Set better y-limits to see the relevant region
ax1_right.set_ylim([0, 8])  # Focus on 0-8 eV range

ax1.set_xlabel('Position x (Å)', fontsize=12)
ax1.set_ylabel('Wavefunction $\psi(x)$', fontsize=12)
ax1.set_title('Quantum Harmonic Oscillator: Wavefunctions with Energy Levels', fontsize=14)
ax1.legend(loc='upper left', fontsize=10)
ax1_right.legend(loc='upper right', fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim([-15, 15])  # Zoom to relevant region

# Plot 2: Probability Densities with potential and energy levels
ax2.plot(x, prob1, label='$|\psi_0|^2$ (n=0)', linewidth=2)
ax2.plot(x, prob2, label='$|\psi_1|^2$ (n=1)', linewidth=2)
ax2.plot(x, prob3, label='$|\psi_2|^2$ (n=2)', linewidth=2)
ax2.plot(x, prob4, label='$|\psi_3|^2$ (n=3)', linewidth=2)
ax2.plot(x, prob5, label='$|\psi_4|^2$ (n=4)', linewidth=2)

# Add potential on secondary y-axis
ax2_right = ax2.twinx()
ax2_right.plot(x, V, 'k--', linewidth=2.5, alpha=0.6, label='V(x)')
ax2_right.set_ylabel('Potential V(x) (eV)', fontsize=12)

# Add horizontal lines for energy levels
for i, (energy, color) in enumerate(zip(E, colors)):
    ax2_right.axhline(y=energy, color=color, linestyle=':', linewidth=1.5, alpha=0.7)

# Set better y-limits
ax2_right.set_ylim([0, 8])

ax2.set_xlabel('Position x (Å)', fontsize=12)
ax2.set_ylabel('Probability Density $|\psi(x)|^2$', fontsize=12)
ax2.set_title('Probability Densities: Peaks at Classical Turning Points', fontsize=14)
ax2.legend(loc='upper left', fontsize=10)
ax2_right.legend(loc='upper right', fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim([-15, 15])

plt.tight_layout()
plt.savefig('results/qho_wavefunctions.png', dpi=300)
plt.show()

print("Plot saved as 'results/qho_wavefunctions.png'")
print("\nNote: Dotted horizontal lines show energy levels E_n.")
print("Probability density peaks occur near where E_n intersects V(x) (classical turning points).")