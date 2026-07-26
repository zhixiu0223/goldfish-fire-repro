"""v5: 牆體構造(材質+厚度)作為可調參數,用文獻 kρc 值比較 石膏板單層/雙層/混凝土"""
import numpy as np
from scipy.integrate import solve_ivp

O2_frac0 = 0.21
O2_mass_frac0 = 0.233
rho_air = 1.2
cp_air = 1.005
O2_heat_release = 13100
heat_of_comb = 16000
c_water = 4.186
T_ambient = 20.0
extinct_hi = 0.21
extinct_lo = 0.12
alpha_fire = 0.1876
DO_lethal = 2.0
T_lethal = 35.0
Q10 = 2.5
base_O2_consumption = 0.05
DO_half_sat = 0.5

# 文獻材料性質 [Drysdale 1999]
WALL_MATERIALS = {
    "gypsum_1layer": dict(k=0.48, rho=1440, c=840, thickness=0.015),   # 單層石膏板 ~15mm
    "gypsum_2layer": dict(k=0.48, rho=1440, c=840, thickness=0.030),   # 雙層石膏板 ~30mm
    "concrete":      dict(k=1.1,  rho=2100, c=880, thickness=0.150),  # 混凝土牆/樓板 ~150mm(k,rho,c取文獻範圍中點)
}

def DO_sat(T):
    T = np.clip(T, 0, 45)
    return 14.6 - 0.41*T + 0.005*T**2

def simulate_full(V_room_m3, R_distance_m, wall_key="gypsum_1layer",
                   V_tank_L=80.0, m_fuel_kg=80.0, t_max=14400,
                   chi_r=0.30, h_conv=20.0, kLa_mult=1.0):
    mat = WALL_MATERIALS[wall_key]
    krhoc = mat["k"]*mat["rho"]*mat["c"]
    alpha_mat = mat["k"]/(mat["rho"]*mat["c"])
    t_valid = mat["thickness"]**2 / (4*alpha_mat)

    V_tank_m3 = V_tank_L / 1000
    side = V_tank_m3 ** (1/3)
    A_tank_total = 6 * side**2
    A_tank_facing = side**2
    m_air = V_room_m3 * rho_air
    A_wall = 6 * (V_room_m3)**(2/3)
    kLa = 0.02 * (A_tank_total / V_tank_m3) / 0.5 * kLa_mult

    def rhs(t, y):
        T_room, O2_frac, fuel_burned, T_water, DO = y
        HRR_growth = alpha_fire * t**2
        fuel_factor = 1.0 if fuel_burned < m_fuel_kg else 0.0
        O2_factor = np.clip((O2_frac - extinct_lo) / (extinct_hi - extinct_lo), 0, 1)
        HRR = min(HRR_growth * fuel_factor * O2_factor, 8000.0)

        O2_mass_rate = HRR / O2_heat_release
        dO2_frac_dt = -O2_mass_rate / (m_air * O2_mass_frac0 / O2_frac0) if m_air > 0 else 0
        dfuel_dt = HRR / heat_of_comb

        t_eff = min(max(t, 1.0), t_valid)
        h_k = np.sqrt(krhoc / (np.pi * t_eff))
        q_room_to_wall = h_k/1000 * A_wall * (T_room - T_ambient)

        dT_room_dt = ((1-chi_r)*HRR - q_room_to_wall) / (m_air * cp_air) if m_air*cp_air > 0 else 0

        q_conv = h_conv/1000 * A_tank_total * (T_room - T_water)
        q_rad_flux = min(chi_r * HRR / (4*np.pi*max(R_distance_m,0.3)**2), 150.0)
        q_rad = q_rad_flux * A_tank_facing
        m_water = V_tank_L
        dT_water_dt = (q_conv + q_rad) / (m_water * c_water)

        DO_eq = DO_sat(T_water) * (O2_frac / O2_frac0)
        d_DO_diffusion = kLa/60 * (DO_eq - DO)
        fish_consumption = (base_O2_consumption/60) * Q10**((T_water-20)/10) * (max(DO,0)/(max(DO,0)+DO_half_sat))
        dDO_dt = d_DO_diffusion - fish_consumption

        return [dT_room_dt, dO2_frac_dt, dfuel_dt, dT_water_dt, dDO_dt]

    y0 = [T_ambient, O2_frac0, 0.0, T_ambient, DO_sat(T_ambient)]
    sol = solve_ivp(rhs, [0, t_max], y0, method='RK45', max_step=5.0)

    T_water_arr = sol.y[3]
    DO_arr = sol.y[4]
    survive = (T_water_arr.max() < T_lethal) and (DO_arr.min() > DO_lethal)
    return survive, T_water_arr.max(), DO_arr.min()

def find_R_crit(V_room, wall_key, R_lo=0.3, R_hi=6.0, tol=0.05):
    s_lo, _, _ = simulate_full(V_room, R_lo, wall_key=wall_key)
    s_hi, _, _ = simulate_full(V_room, R_hi, wall_key=wall_key)
    if s_lo and s_hi: return R_lo
    if not s_lo and not s_hi: return R_hi
    while R_hi - R_lo > tol:
        R_mid = (R_lo+R_hi)/2
        s_mid, _, _ = simulate_full(V_room, R_mid, wall_key=wall_key)
        if s_mid: R_hi = R_mid
        else: R_lo = R_mid
    return (R_lo+R_hi)/2

if __name__ == "__main__":
    for key, mat in WALL_MATERIALS.items():
        krhoc = mat["k"]*mat["rho"]*mat["c"]
        alpha_mat = mat["k"]/(mat["rho"]*mat["c"])
        t_valid = mat["thickness"]**2/(4*alpha_mat)
        print(f"{key:16s} kρc={krhoc:.3e}  厚度={mat['thickness']*1000:.0f}mm  半無限假設有效時間={t_valid/60:.1f} 分鐘")
