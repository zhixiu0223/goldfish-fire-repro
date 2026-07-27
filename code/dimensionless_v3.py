"""
無因次分析v3:修正t_char——原本是固定常數,沒有隨房間體積(氧氣供應量)縮放,
導致大房間案例的輻射曝露時間被低估,是NR-04近距離誤分類的真正根因。
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dimensionless import simulate_with_pi
from dimensionless_v2 import max_fuel_burnable_by_O2
import fish_ode_model5 as m5
from fish_ode_model5 import DO_sat, T_ambient, T_lethal, DO_lethal

def compute_Pi_v3(V_room_m3, R_distance_m, V_tank_L=80.0, m_fuel_kg=80.0,
                    chi_r=0.30, h_conv=20.0):
    m_burn_O2limited = max_fuel_burnable_by_O2(V_room_m3)
    m_burn = min(m_burn_O2limited, m_fuel_kg)
    E_fire = m_burn * m5.heat_of_comb  # kJ

    # 修正: t_char改用「以t^2成長曲線燒完E_fire總能量所需時間」,隨房間體積(氧氣量)縮放
    # E_fire = integral(alpha*t^2)dt from 0 to t_burn = alpha*t_burn^3/3
    t_char = (3*E_fire / m5.alpha_fire) ** (1/3) if E_fire > 0 else 1.0

    V_tank_m3 = V_tank_L / 1000
    side = V_tank_m3 ** (1/3)
    A_facing = side**2

    q_rad_flux_char = chi_r * (E_fire/t_char) / (4*np.pi*max(R_distance_m,0.3)**2)
    q_rad_flux_char = min(q_rad_flux_char, 150.0)  # 跟實際ODE一致的封頂
    Q_rad_est = min(q_rad_flux_char * A_facing * t_char, chi_r*E_fire)

    Q_conv_est = (1-chi_r) * E_fire * 0.3  # 仍是未校準的0.3係數,見Future Work
    Q_total_est = Q_rad_est + Q_conv_est

    Pi_heat = Q_total_est / (V_tank_L * m5.c_water * (T_lethal - T_ambient))
    return Pi_heat

V_room_range = np.linspace(20, 300, 12)
R_dist_range = np.linspace(0.3, 4.0, 12)

results = []
for V in V_room_range:
    for R in R_dist_range:
        r = simulate_with_pi(V, R)
        pi_v3 = compute_Pi_v3(V, R)
        results.append((V, R, r["survive"], pi_v3))

results = np.array(results)
survive = results[:,2].astype(bool)
Pi = results[:,3]

thresholds = np.linspace(Pi.min(), Pi.max(), 300)
best_acc, best_thresh = 0, 0
for th in thresholds:
    acc = np.mean((Pi >= th) == ~survive)
    if acc > best_acc:
        best_acc, best_thresh = acc, th

print(f"v3(t_char隨V_room縮放)準確率: {best_acc*100:.1f}%, 門檻={best_thresh:.3f}")
mis = (Pi >= best_thresh) != ~survive
print(f"誤分類數: {mis.sum()}")
if mis.sum() > 0:
    print("誤分類樣本:")
    for i in np.where(mis)[0]:
        print(f"  V={results[i,0]:.1f} R={results[i,1]:.2f} 實際={'死' if not survive[i] else '活'} Pi={Pi[i]:.2f}")
