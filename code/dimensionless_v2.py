"""
修正版無因次分析:Pi群必須完全由「跑模擬前就知道的輸入參數」構成,
不能用模擬輸出(如實際水溫、實際累積熱量)去反推,否則是循環論證。
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dimensionless import simulate_with_pi
import fish_ode_model5 as m5
from fish_ode_model5 import DO_sat, T_ambient, T_lethal, DO_lethal

def max_fuel_burnable_by_O2(V_room_m3, extinct_frac=m5.extinct_lo):
    """純用輸入參數(房間體積)算:氧氣濃度掉到熄滅門檻前,理論上最多燒掉多少燃料"""
    n_air = V_room_m3 * 1000 / 22.4
    n_O2_avail = n_air * (m5.O2_frac0 - extinct_frac)
    mass_O2 = n_O2_avail * 32 / 1000
    wood_O2_demand = 1.2
    return mass_O2 / wood_O2_demand

def compute_Pi_input_only(V_room_m3, R_distance_m, V_tank_L=80.0, m_fuel_kg=80.0,
                            chi_r=0.30, h_conv=20.0):
    """完全用輸入參數建構的無因次群,不碰任何模擬輸出"""
    # 實際會燒掉的燃料量:氧氣限制 vs 燃料總量,取較小值(純輸入參數推導,不跑ODE)
    m_burn_O2limited = max_fuel_burnable_by_O2(V_room_m3)
    m_burn = min(m_burn_O2limited, m_fuel_kg)
    E_fire = m_burn * m5.heat_of_comb  # kJ,理論總放熱量上限

    V_tank_m3 = V_tank_L / 1000
    side = V_tank_m3 ** (1/3)
    A_facing = side**2
    A_tank_total = 6*side**2

    # 輻射分量(有R依賴):輻射熱通量估計乘上一個粗略的"有效曝露時間"尺度
    # 用火災成長曲線的特徵時間(HRR達到有意義量級的時間尺度)去估計曝露時間,而非模擬結果
    # t_char: alpha*t^2 = 某特徵功率(例如 1000kW)時的時間,當作曝露時間尺度
    t_char = np.sqrt(1000.0 / m5.alpha_fire)  # 秒

    Q_rad_est = chi_r * E_fire * A_facing / (4*np.pi*max(R_distance_m,0.3)**2) / (E_fire/1) * min(E_fire, chi_r*E_fire)  # 先簡化見下

    # 簡化重新定義:輻射能量估計 = 輻射通量(用E_fire及R估計的特徵通量) x 特徵曝露時間,上限為E_fire本身
    q_rad_flux_char = chi_r * (E_fire/t_char) / (4*np.pi*max(R_distance_m,0.3)**2)  # kW/m2 概估
    Q_rad_est = min(q_rad_flux_char * A_facing * t_char, chi_r*E_fire)  # kJ,不超過輻射能量上限

    # 對流分量(0D模型下不隨R變化,只隨房間規模/燃料量變化)
    Q_conv_est = (1-chi_r) * E_fire * 0.3  # 概估:非輻射熱量中,有多少比例能傳到水缸(拍腦袋係數0.3,需標註不確定)

    Q_total_est = Q_rad_est + Q_conv_est

    Pi_heat_apriori = Q_total_est / (V_tank_L * m5.c_water * (T_lethal - T_ambient))
    return Pi_heat_apriori

V_room_range = np.linspace(20, 300, 12)
R_dist_range = np.linspace(0.3, 4.0, 12)

results = []
for V in V_room_range:
    for R in R_dist_range:
        r = simulate_with_pi(V, R)  # 真實跑模擬拿到survive/death(ground truth)
        pi_apriori = compute_Pi_input_only(V, R)  # 純輸入參數算的Pi群(預測用)
        results.append((V, R, r["survive"], pi_apriori))

results_arr = np.array([(v,r,int(s),p) for v,r,s,p in results])
survive = results_arr[:,2].astype(bool)
Pi_apriori = results_arr[:,3]

# 找最佳分類門檻
thresholds = np.linspace(Pi_apriori.min(), Pi_apriori.max(), 200)
best_acc, best_thresh = 0, 0
for th in thresholds:
    pred_death = Pi_apriori >= th
    acc = np.mean(pred_death == ~survive)
    if acc > best_acc:
        best_acc, best_thresh = acc, th

print(f"死亡案例 Pi_apriori 範圍: {Pi_apriori[~survive].min():.4f} ~ {Pi_apriori[~survive].max():.4f}")
print(f"存活案例 Pi_apriori 範圍: {Pi_apriori[survive].min():.4f} ~ {Pi_apriori[survive].max():.4f}")
print(f"最佳單一門檻值 = {best_thresh:.4f}, 分類準確率 = {best_acc*100:.1f}%")

fig, ax = plt.subplots(figsize=(9,4))
ax.scatter(Pi_apriori[survive], np.zeros(np.sum(survive))+0.05, c='green', label='Survive (simulated)', alpha=0.6, s=30)
ax.scatter(Pi_apriori[~survive], np.zeros(np.sum(~survive))-0.05, c='red', label='Death (simulated)', alpha=0.6, s=30)
ax.axvline(best_thresh, color='gray', ls='--', label=f'best threshold={best_thresh:.3f}')
ax.set_xlabel('Pi_heat (a priori, from input parameters only)')
ax.set_ylim(-0.3,0.3); ax.set_yticks([])
ax.set_title(f'Does an INPUT-ONLY Pi group predict survive/death? (accuracy={best_acc*100:.0f}%)')
ax.legend()
plt.tight_layout()
plt.savefig("dimensionless_v2.png", dpi=150)
print("saved")
