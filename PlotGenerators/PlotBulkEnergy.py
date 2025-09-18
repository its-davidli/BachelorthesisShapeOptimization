import numpy as np
import matplotlib.pyplot as plt

# Parameters (choose values as needed)
U0 = 0
B = -1 # More negative B increases the cubic term's effect
C = 1    # Smaller C sharpens the quartic term

A_Star = (B**2)/(24*C)
A_c = (B**2)/(27*C)


S = np.linspace(-0.15, 0.5, 400)
# psi_bulk = U0 + (A/3)*S**2 + (2*B/27)*S**3 + (C/9)*S**4

plt.figure(figsize=(6,4))

A_values = [A_c/1.1, A_c ,A_c + (A_Star- A_c)/2]
labels = [r'$A<A_c$', r'$A=A_c$', r'$A>A_c$']
for A, label in zip(A_values, labels):
    psi_bulk = U0 + (A/3)*S**2 + (2*B/27)*S**3 + (C/9)*S**4
    plt.plot(S, psi_bulk, label=fr'$\psi_\mathrm{{bulk}}(S)$, {label}')

    # Mark the point (-B + sqrt(B^2 - 24AC)) / (4C) for each A
    for A in A_values:
        discriminant = B**2 - 24*A*C
        if discriminant >= 0:
            S_crit = (-B + np.sqrt(discriminant)) / (4*C)
            psi_crit = U0 + (A/3)*S_crit**2 + (2*B/27)*S_crit**3 + (C/9)*S_crit**4
            plt.plot(S_crit, psi_crit, 'o', markersize=3)
plt.xlabel(r'$S$')
plt.ylabel(r'Bulk energy density $\psi_\mathrm{bulk}$')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.legend()
plt.tight_layout()
plt.yticks([0], [r'$0$'])
plt.show()