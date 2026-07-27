"""
追NR-04的線索:誤分類樣本全集中在R<=0.64m,是不是因為Pi_heat公式裡的
輻射估計項,跟實際ODE模擬裡真正發生的輻射通量算法不一致(尤其是近距離的封頂處理)?
"""
import numpy as np
from dimensionless import simulate_with_pi
from dimensionless_v2 import compute_Pi_input_only
import fish_ode_model5 as m5

# 針對誤分類的其中一個代表案例,實際比對
V_test, R_test = 147.3, 0.30

# a priori公式估計的輻射能量
def diagnose(V, R, chi_r=0.30):
    E_fire = min(m5.heat_of_comb * 80.0,  # m_fuel_kg上限
                 m5.heat_of_comb * (lambda v: min(v * 1000/22.4*(0.21-0.12)*32/1000/1.2, 80.0))(V))
    t_char = np.sqrt(1000.0 / m5.alpha_fire)
    q_rad_flux_char = chi_r * (E_fire/t_char) / (4*np.pi*max(R,0.3)**2)
    print(f"V={V}, R={R}")
    print(f"  a priori 特徵輻射通量 q_rad_flux_char = {q_rad_flux_char:.1f} kW/m^2")
    print(f"  實際ODE模擬裡的輻射通量上限(硬編碼) = 150.0 kW/m^2")
    if q_rad_flux_char > 150:
        print(f"  >>> a priori公式沒有套用150上限,估計值被嚴重高估 {q_rad_flux_char/150:.1f}倍")
    else:
        print(f"  未觸及上限")

diagnose(147.3, 0.30)
print()
diagnose(96.4, 0.30)
print()
diagnose(274.5, 0.64)
