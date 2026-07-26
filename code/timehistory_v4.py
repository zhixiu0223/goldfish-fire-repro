import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fish_ode_model4 import simulate_full, DO_lethal, T_lethal

cases = {
    "Clear survive\n(V=50 m3, R=3.0m)": simulate_full(50, 3.0, t_max=14400),
    "Near boundary\n(V=150 m3, R=0.9m)": simulate_full(150, 0.9, t_max=14400),
    "Clear death\n(V=250 m3, R=0.5m)": simulate_full(250, 0.5, t_max=14400),
}

fig, axes = plt.subplots(3, 3, figsize=(16, 9), sharex='col')
for col, (title, r) in enumerate(cases.items()):
    t_hr = r["t"] / 3600
    axes[0, col].plot(t_hr, r["T_room"], color='orange')
    axes[0, col].set_title(title, fontsize=11); axes[0, col].set_ylabel("T_room (C)")
    axes[1, col].plot(t_hr, r["T_water"], color='crimson')
    axes[1, col].axhline(T_lethal, color='red', ls='--', lw=1, label=f'lethal {T_lethal}C')
    axes[1, col].set_ylabel("T_water (C)"); axes[1, col].legend(fontsize=8)
    axes[2, col].plot(t_hr, r["DO"], color='green')
    axes[2, col].axhline(DO_lethal, color='red', ls='--', lw=1, label=f'lethal {DO_lethal}mg/L')
    axes[2, col].set_ylabel("DO (mg/L)"); axes[2, col].set_xlabel("Time (hr)")
    axes[2, col].legend(fontsize=8)

plt.suptitle("v4: literature-calibrated wall heat loss (semi-infinite gypsum conduction)", fontsize=12)
plt.tight_layout()
plt.savefig("timehistory_v4.png", dpi=140)
print("saved")

for title, r in cases.items():
    print(f"\n=== {title.splitlines()[0]} ===")
    print(f"  final T_room={r['T_room'][-1]:.1f}C  T_water={r['T_water'][-1]:.1f}C  DO={r['DO'][-1]:.2f}mg/L")
    print(f"  max T_water={r['T_water'].max():.1f}C at t={r['t'][np.argmax(r['T_water'])]/60:.1f}min")
    print(f"  min DO={r['DO'].min():.2f}mg/L at t={r['t'][np.argmin(r['DO'])]/60:.1f}min")
    print(f"  crosses heat-lethal? {(r['T_water']>=T_lethal).any()}   crosses DO-lethal? {(r['DO']<=DO_lethal).any()}")
