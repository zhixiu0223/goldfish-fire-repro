import numpy as np
from dimensionless import simulate_with_pi
from dimensionless_v2 import compute_Pi_input_only

V_room_range = np.linspace(20, 300, 12)
R_dist_range = np.linspace(0.3, 4.0, 12)

results = []
for V in V_room_range:
    for R in R_dist_range:
        r = simulate_with_pi(V, R)
        pi = compute_Pi_input_only(V, R)
        results.append((V, R, r["survive"], pi))

results = np.array(results, dtype=float)
survive = results[:,2].astype(bool)
Pi = results[:,3]

best_thresh = 22.0
pred_death = Pi >= best_thresh
actual_death = ~survive
misclassified = pred_death != actual_death

print(f"總樣本數: {len(results)}, 誤分類數: {misclassified.sum()}")
print("\n誤分類的樣本:")
print(f"{'V_room':>8} {'R':>6} {'實際':>6} {'預測':>6} {'Pi_heat':>8}")
for i in np.where(misclassified)[0]:
    V, R, s, p = results[i]
    actual = "存活" if s else "死亡"
    pred = "死亡" if pred_death[i] else "存活"
    print(f"{V:8.1f} {R:6.2f} {actual:>6} {pred:>6} {p:8.3f}")

print(f"\n誤分類樣本的 V_room 範圍: {results[misclassified,0].min():.1f} ~ {results[misclassified,0].max():.1f}")
print(f"誤分類樣本的 R 範圍: {results[misclassified,1].min():.2f} ~ {results[misclassified,1].max():.2f}")
print(f"誤分類樣本的 Pi_heat 範圍: {Pi[misclassified].min():.3f} ~ {Pi[misclassified].max():.3f}")
