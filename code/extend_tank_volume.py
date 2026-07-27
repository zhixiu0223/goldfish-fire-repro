"""
把水缸體積也納入同一組collapse測試,看Pi_heat(v2,95.1%準確率版本)能不能同時處理這個新維度
"""
import numpy as np
from dimensionless import simulate_with_pi
from dimensionless_v2 import compute_Pi_input_only

V_room_range = np.linspace(30, 250, 8)
R_dist_range = np.linspace(0.3, 3.0, 8)
V_tank_range = [40.0, 80.0, 160.0]  # 小/中/大水缸

results = []
for V_tank in V_tank_range:
    for V in V_room_range:
        for R in R_dist_range:
            r = simulate_with_pi(V, R, V_tank_L=V_tank)
            pi = compute_Pi_input_only(V, R, V_tank_L=V_tank)
            results.append((V, R, V_tank, r["survive"], pi))

results = np.array(results)
survive = results[:,3].astype(bool)
Pi = results[:,4]

thresholds = np.linspace(Pi.min(), Pi.max(), 300)
best_acc, best_thresh = 0, 0
for th in thresholds:
    acc = np.mean((Pi >= th) == ~survive)
    if acc > best_acc:
        best_acc, best_thresh = acc, th

print(f"加入水缸體積(40/80/160L)後,同一個Pi_heat公式分類準確率: {best_acc*100:.1f}%, 門檻={best_thresh:.3f}")
print(f"總樣本數: {len(results)}")

# 分開看每個水缸體積的準確率,看水缸體積是不是新的破口
for V_tank in V_tank_range:
    mask = results[:,2] == V_tank
    sub_survive = survive[mask]
    sub_Pi = Pi[mask]
    pred = sub_Pi >= best_thresh
    acc = np.mean(pred == ~sub_survive)
    print(f"  水缸{V_tank}L 子集準確率: {acc*100:.1f}% (n={mask.sum()})")
