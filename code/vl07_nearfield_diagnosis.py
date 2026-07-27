"""
VL-07: 把實際ODE模擬裡的Q_delivered拆成輻射/對流兩個分量,
跟a priori公式裡的Q_rad_est/Q_conv_est分別比對,
找出近距離(R<=0.64m)誤分類究竟是哪一項被低估。
"""
import numpy as np
from scipy.integrate import solve_ivp
import fish_ode_model5 as m5
from fish_ode_model5 import DO_sat, T_ambient, T_lethal, DO_lethal
from dimensionless_v2 import max_fuel_burnable_by_O2

def simulate_decomposed(V_room_m3, R_distance_m, wall_key="gypsum_1layer",
                         V_tank_L=80.0, m_fuel_kg=80.0, t_max=14400,
                         chi_r=0.30, h_conv=20.0, kLa_mult=1.0):
    mat = m5.WALL_MATERIALS[wall_key]
    krhoc = mat["k"]*mat["rho"]*mat["c"]
    alpha_mat = mat["k"]/(mat["rho"]*mat["c"])
    t_valid = mat["thickness"]**2 / (4*alpha_mat)

    V_tank_m3 = V_tank_L / 1000
    side = V_tank_m3 ** (1/3)
    A_tank_total = 6 * side**2
    A_tank_facing = side**2
    m_air = V_room_m3 * m5.rho_air
    A_wall = 6 * (V_room_m3)**(2/3)
    kLa = 0.02 * (A_tank_total / V_tank_m3) / 0.5 * kLa_mult

    def rhs(t, y):
        T_room, O2_frac, fuel_burned, T_water, DO, Q_rad_cum, Q_conv_cum = y
        HRR_growth = m5.alpha_fire * t**2
        fuel_factor = 1.0 if fuel_burned < m_fuel_kg else 0.0
        O2_factor = np.clip((O2_frac - m5.extinct_lo) / (m5.extinct_hi - m5.extinct_lo), 0, 1)
        HRR = min(HRR_growth * fuel_factor * O2_factor, 8000.0)

        O2_mass_rate = HRR / m5.O2_heat_release
        dO2_frac_dt = -O2_mass_rate / (m_air * m5.O2_mass_frac0 / m5.O2_frac0) if m_air > 0 else 0
        dfuel_dt = HRR / m5.heat_of_comb

        t_eff = min(max(t, 1.0), t_valid)
        h_k = np.sqrt(krhoc / (np.pi * t_eff))
        q_room_to_wall = h_k/1000 * A_wall * (T_room - T_ambient)
        dT_room_dt = ((1-chi_r)*HRR - q_room_to_wall) / (m_air * m5.cp_air) if m_air*m5.cp_air > 0 else 0

        q_conv = h_conv/1000 * A_tank_total * (T_room - T_water)
        q_rad_flux = min(chi_r * HRR / (4*np.pi*max(R_distance_m,0.3)**2), 150.0)
        q_rad = q_rad_flux * A_tank_facing
        m_water = V_tank_L
        dT_water_dt = (q_conv + q_rad) / (m_water * m5.c_water)

        DO_eq = DO_sat(T_water) * (O2_frac / m5.O2_frac0)
        d_DO_diffusion = kLa/60 * (DO_eq - DO)
        fish_consumption = (m5.base_O2_consumption/60) * m5.Q10**((T_water-20)/10) * (max(DO,0)/(max(DO,0)+m5.DO_half_sat))
        dDO_dt = d_DO_diffusion - fish_consumption

        return [dT_room_dt, dO2_frac_dt, dfuel_dt, dT_water_dt, dDO_dt, q_rad, q_conv]

    y0 = [T_ambient, m5.O2_frac0, 0.0, T_ambient, DO_sat(T_ambient), 0.0, 0.0]
    sol = solve_ivp(rhs, [0, t_max], y0, method='RK45', max_step=5.0)

    T_water_arr = sol.y[3]
    DO_arr = sol.y[4]
    Q_rad_actual = sol.y[5][-1]
    Q_conv_actual = sol.y[6][-1]
    survive = (T_water_arr.max() < T_lethal) and (DO_arr.min() > DO_lethal)
    return survive, Q_rad_actual, Q_conv_actual

def compute_Pi_components(V_room_m3, R_distance_m, V_tank_L=80.0, m_fuel_kg=80.0, chi_r=0.30):
    m_burn_O2limited = max_fuel_burnable_by_O2(V_room_m3)
    m_burn = min(m_burn_O2limited, m_fuel_kg)
    E_fire = m_burn * m5.heat_of_comb
    t_char = np.sqrt(1000.0 / m5.alpha_fire)
    V_tank_m3 = V_tank_L / 1000
    side = V_tank_m3 ** (1/3)
    A_facing = side**2
    q_rad_flux_char = chi_r * (E_fire/t_char) / (4*np.pi*max(R_distance_m,0.3)**2)
    Q_rad_est = min(q_rad_flux_char * A_facing * t_char, chi_r*E_fire)
    Q_conv_est = (1-chi_r) * E_fire * 0.3
    return Q_rad_est, Q_conv_est

print(f"{'V_room':>7} {'R':>5} {'實際':>5} {'Q_rad實際':>10} {'Q_rad估計':>10} {'比值':>6} | {'Q_conv實際':>11} {'Q_conv估計':>11} {'比值':>6}")
misclassified_cases = [(96.4,0.30),(147.3,0.30),(198.2,0.30),(274.5,0.30),(274.5,0.64)]
for V, R in misclassified_cases:
    survive, Qr_act, Qc_act = simulate_decomposed(V, R)
    Qr_est, Qc_est = compute_Pi_components(V, R)
    actual = "死" if not survive else "活"
    ratio_r = Qr_act/Qr_est if Qr_est>0 else float('inf')
    ratio_c = Qc_act/Qc_est if Qc_est>0 else float('inf')
    print(f"{V:7.1f} {R:5.2f} {actual:>5} {Qr_act:10.1f} {Qr_est:10.1f} {ratio_r:6.2f} | {Qc_act:11.1f} {Qc_est:11.1f} {ratio_c:6.2f}")
