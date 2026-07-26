import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fish_ode_model5 import WALL_MATERIALS, DO_sat, T_ambient, T_lethal, DO_lethal
import fish_ode_model5 as m5

def simulate_with_pi(V_room_m3, R_distance_m, wall_key="gypsum_1layer",
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

    Q_delivered_holder = {"Q": 0.0}  # kJ, 用簡單trick在rhs外部累積

    def rhs(t, y):
        T_room, O2_frac, fuel_burned, T_water, DO, Q_cum = y
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
        dQ_cum_dt = q_conv + q_rad  # kW -> 累積 kJ

        DO_eq = DO_sat(T_water) * (O2_frac / m5.O2_frac0)
        d_DO_diffusion = kLa/60 * (DO_eq - DO)
        fish_consumption = (m5.base_O2_consumption/60) * m5.Q10**((T_water-20)/10) * (max(DO,0)/(max(DO,0)+m5.DO_half_sat))
        dDO_dt = d_DO_diffusion - fish_consumption

        return [dT_room_dt, dO2_frac_dt, dfuel_dt, dT_water_dt, dDO_dt, dQ_cum_dt]

    y0 = [T_ambient, m5.O2_frac0, 0.0, T_ambient, DO_sat(T_ambient), 0.0]
    sol = solve_ivp(rhs, [0, t_max], y0, method='RK45', max_step=5.0)

    T_water_arr = sol.y[3]
    DO_arr = sol.y[4]
    Q_delivered = sol.y[5][-1]  # kJ, 累積傳給水的總熱量

    max_T = T_water_arr.max()
    min_DO = DO_arr.min()
    survive = (max_T < T_lethal) and (min_DO > DO_lethal)

    m_water = V_tank_L
    Pi_heat = Q_delivered / (m_water * m5.c_water * (T_lethal - T_ambient))
    DO0 = DO_sat(T_ambient)
    Pi_hypoxia = (DO0 - min_DO) / (DO0 - DO_lethal)

    return dict(survive=survive, max_T=max_T, min_DO=min_DO,
                Pi_heat=Pi_heat, Pi_hypoxia=Pi_hypoxia)

if __name__ == "__main__":
    r = simulate_with_pi(150, 0.9)
    print(r)
