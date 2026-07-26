import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fish_ode_model5 import find_R_crit, WALL_MATERIALS

V_room_range = np.linspace(20, 300, 12)

fig, ax = plt.subplots(figsize=(9, 6))
colors = {"gypsum_1layer": "tab:blue", "gypsum_2layer": "tab:orange", "concrete": "tab:green"}
labels = {"gypsum_1layer": "Single-layer gypsum (15mm)",
          "gypsum_2layer": "Double-layer gypsum (30mm)",
          "concrete": "Concrete (150mm)"}

for wall_key in WALL_MATERIALS:
    R_crits = [find_R_crit(V, wall_key) for V in V_room_range]
    ax.plot(V_room_range, R_crits, marker='o', color=colors[wall_key], label=labels[wall_key])
    print(f"{wall_key:16s} R_crit: {[f'{r:.2f}' for r in R_crits]}")

ax.set_xlabel("Room volume (m^3)", fontsize=12)
ax.set_ylabel("Critical distance R_crit (m)\n(above = survive, below = death)", fontsize=12)
ax.set_title("Sensitivity of survival boundary to wall construction", fontsize=13)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("wall_sensitivity.png", dpi=150)
print("saved")
