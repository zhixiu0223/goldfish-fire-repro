"""
系統性掃描near-field失效邊界:對每個V_room,細掃R從0.3到1.5m,
找出Pi_heat(門檻22.0)開始跟實際模擬一致的臨界距離 R_fail(V_room)。
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dimensionless import simulate_with_pi
from dimensionless_v2 import compute_Pi_input_only

V_room_range = np.linspace(20, 300, 15)
R_fine = np.linspace(0.3, 1.5, 25)  # 密集掃描near-field區域
THRESH = 22.0

results = []
for V in V_room_range:
    for R in R_fine:
        r = simulate_with_pi(V, R)
        pi = compute_Pi_input_only(V, R)
        pred_death = pi >= THRESH
        actual_death = not r["survive"]
        correct = pred_death == actual_death
        results.append((V, R, actual_death, pred_death, correct))

results = np.array(results)

# 對每個V_room,找出「從R小到大,連續正確分類開始的最小R」
R_fail_boundary = []
for V in V_room_range:
    mask = results[:,0] == V
    sub = results[mask]
    sub = sub[np.argsort(sub[:,1])]  # 按R排序
    incorrect_Rs = sub[sub[:,4]==0, 1]
    if len(incorrect_Rs) > 0:
        R_fail_boundary.append((V, incorrect_Rs.max()))  # 最大的錯誤發生點
    else:
        R_fail_boundary.append((V, None))  # 這個V_room完全沒有錯誤

print(f"{'V_room':>8} {'最大誤判R':>10} {'該V_room誤判樣本數':>18}")
total_wrong = 0
for V, R_max_wrong in R_fail_boundary:
    mask = results[:,0]==V
    n_wrong = int(np.sum(results[mask,4]==0))
    total_wrong += n_wrong
    print(f"{V:8.1f} {str(R_max_wrong) if R_max_wrong else '無誤判':>10} {n_wrong:18d}")

print(f"\n總誤判數: {total_wrong} / {len(results)} = {total_wrong/len(results)*100:.1f}%")

# 畫圖:每個V_room的誤判範圍
fig, ax = plt.subplots(figsize=(9,6))
for V in V_room_range:
    mask = results[:,0]==V
    sub = results[mask]
    correct_R = sub[sub[:,4]==1, 1]
    wrong_R = sub[sub[:,4]==0, 1]
    ax.scatter([V]*len(correct_R), correct_R, c='green', s=15, alpha=0.5)
    ax.scatter([V]*len(wrong_R), wrong_R, c='red', s=25, marker='x')

ax.scatter([],[], c='green', label='Pi_heat correct')
ax.scatter([],[], c='red', marker='x', label='Pi_heat wrong (misclassified)')
ax.set_xlabel('Room volume (m^3)')
ax.set_ylabel('Distance R (m)')
ax.set_title('Near-field failure boundary: fine-grained scan (R: 0.3-1.5m)')
ax.legend()
plt.tight_layout()
plt.savefig('nearfield_boundary.png', dpi=150)
print("\nsaved")
