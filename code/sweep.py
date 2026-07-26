import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fish_ode_model import simulate

V_room_range = np.linspace(15, 300, 14)
R_dist_range = np.linspace(0.3, 4.0, 14)

cause_map = np.zeros((len(R_dist_range), len(V_room_range)), dtype=int)
# 0=存活, 1=缺氧死, 2=熱死, 3=熱+缺氧接近同時

cause_to_code = {"存活": 0, "缺氧死": 1, "熱死": 2, "熱+缺氧(接近同時)": 3}

results = []
for i, R in enumerate(R_dist_range):
    for j, V in enumerate(V_room_range):
        r = simulate(V_room_m3=V, R_distance_m=R, V_tank_L=80.0, m_fuel_kg=80.0, t_max=2400)
        cause_map[i, j] = cause_to_code[r["cause"]]
        results.append((V, R, r["cause"], r["max_T_water"], r["min_DO"]))
        print(f"V_room={V:6.1f} R={R:.2f}  -> {r['cause']:20s} maxT={r['max_T_water']:.1f} minDO={r['min_DO']:.2f}")

np.save("cause_map.npy", cause_map)
np.save("V_room_range.npy", V_room_range)
np.save("R_dist_range.npy", R_dist_range)

# 畫相圖
fig, ax = plt.subplots(figsize=(10, 7))
cmap = matplotlib.colors.ListedColormap(['#4CAF50', '#FFC107', '#F44336', '#9C27B0'])
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

im = ax.imshow(cause_map, origin='lower', aspect='auto', cmap=cmap, norm=norm,
               extent=[V_room_range.min(), V_room_range.max(), R_dist_range.min(), R_dist_range.max()])

ax.set_xlabel('房間體積 V_room (m³)', fontsize=12)
ax.set_ylabel('水缸-火源距離 R (m)', fontsize=12)
ax.set_title('金魚存活相圖 (水缸80L, 家具80kg, "極其猛烈"火災成長曲線)', fontsize=13)

cbar = fig.colorbar(im, ax=ax, ticks=[0,1,2,3])
cbar.ax.set_yticklabels(['存活', '缺氧死', '熱死', '熱+缺氧同時'])

plt.tight_layout()
plt.savefig('phase_diagram.png', dpi=150)
print("\n圖已存檔")
