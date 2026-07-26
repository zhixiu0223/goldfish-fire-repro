"""
Collapse analysis: 檢查 Pi_heat(純輸入參數版)的等高線,
是否跟原本2D相圖(V_room x R)裡真實模擬出來的存活/死亡邊界重合。
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dimensionless import simulate_with_pi
from dimensionless_v2 import compute_Pi_input_only

V_room_range = np.linspace(15, 300, 40)
R_dist_range = np.linspace(0.3, 4.0, 40)

VV, RR = np.meshgrid(V_room_range, R_dist_range)
Pi_grid = np.zeros_like(VV)
for i in range(VV.shape[0]):
    for j in range(VV.shape[1]):
        Pi_grid[i,j] = compute_Pi_input_only(VV[i,j], RR[i,j])

# 真實模擬結果(較粗網格,節省運算)
V_coarse = np.linspace(15, 300, 16)
R_coarse = np.linspace(0.3, 4.0, 16)
sim_points = []
for V in V_coarse:
    for R in R_coarse:
        r = simulate_with_pi(V, R)
        sim_points.append((V, R, r["survive"]))

fig, ax = plt.subplots(figsize=(9,7))

# 背景:Pi_heat 等高線
levels = [5, 10, 15, 20, 22, 25, 30, 40, 60]
cs = ax.contour(VV, RR, Pi_grid, levels=levels, cmap='coolwarm', linewidths=1)
ax.clabel(cs, inline=True, fontsize=8, fmt='Π=%.0f')

# 疊加真實模擬的存活/死亡散點
for V, R, survive in sim_points:
    ax.scatter(V, R, c='green' if survive else 'red', s=60,
               marker='o' if survive else 'x', zorder=5)

ax.scatter([], [], c='green', marker='o', label='Simulated: Survive')
ax.scatter([], [], c='red', marker='x', label='Simulated: Death')
ax.plot([], [], color='gray', label='Π_heat contours (input-only formula)')

ax.set_xlabel('Room volume (m^3)')
ax.set_ylabel('Tank-fire distance (m)')
ax.set_title('Collapse check: does a single Π_heat contour ≈ true survive/death boundary?')
ax.legend(loc='upper right', fontsize=9)
plt.tight_layout()
plt.savefig('collapse_analysis.png', dpi=150)
print("saved")

# 量化:最佳門檻contour跟實際邊界的吻合程度(用之前算過的95.1%)
