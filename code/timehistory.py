import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 複用 fish_ode_model2 的常數與物理式,但這次要保留完整時間序列
from fish_ode_model2 import (O2_frac0, O2_mass_frac0, rho_air, cp_air, O2_heat_release,
                              heat_of_comb, h_wall, wall_thick_eff, rho_wall, cp_wall,
                              c_water, T_ambient, extinct_hi, extinct_lo, alpha_fire,
                              DO_lethal, T_lethal, Q10, base_O2_consumption, DO_sat)

def simulate_full(V_room_m3, R_distance_m, V_tank_L=80.0, m_fuel_kg=80.0, t_max=2400,
                   chi_r=0.30, h_conv=20.0, kLa_mult=1.0):
    V_tank_m3 = V_tank_L / 1000
    side = V_tank_m3 ** (1/3)
    A_tank_total = 6 * side**2
    A_tank_facing = side**2
    m_air = V_room_m3 * rho_air
    A_wall = 6 * (V_room_m3)**(2/3)
    m_wall = A_wall * wall_thick_eff * rho_wall
    kLa = 0.02 * (A_tank_total / V_tank_m3) / 0.5 * kLa_mult

    HRR_history = []

    def rhs(t, y):
        T_room, O2_frac, T_wall, fuel_burned, T_water, DO = y
        HRR_growth = alpha_fire * t**2
        fuel_factor = 1.0 if fuel_burned < m_fuel_kg else 0.0
        O2_factor = np.clip((O2_frac - extinct_lo) / (extinct_hi - extinct_lo), 0, 1)
        HRR = min(HRR_growth * fuel_factor * O2_factor, 8000.0)
        HRR_history.append((t, HRR))

        O2_mass_rate = HRR / O2_heat_release
        dO2_frac_dt = -O2_mass_rate / (m_air * O2_mass_frac0 / O2_frac0) if m_air > 0 else 0
        dfuel_dt = HRR / heat_of_comb

        q_room_to_wall = h_wall/1000 * A_wall * (T_room - T_wall)
        dT_room_dt = ((1-chi_r)*HRR - q_room_to_wall) / (m_air * cp_air) if m_air*cp_air > 0 else 0
        dT_wall_dt = q_room_to_wall / (m_wall * cp_wall) if m_wall*cp_wall > 0 else 0

        q_conv = h_conv/1000 * A_tank_total * max(T_room - T_water, 0)
        q_rad_flux = min(chi_r * HRR / (4*np.pi*max(R_distance_m,0.3)**2), 150.0)
        q_rad = q_rad_flux * A_tank_facing
        m_water = V_tank_L
        dT_water_dt = (q_conv + q_rad) / (m_water * c_water)

        DO_eq = DO_sat(T_water) * (O2_frac / O2_frac0)
        d_DO_diffusion = kLa/60 * (DO_eq - DO)
        fish_consumption = (base_O2_consumption/60) * Q10**((T_water-20)/10)
        dDO_dt = d_DO_diffusion - fish_consumption

        return [dT_room_dt, dO2_frac_dt, dT_wall_dt, dfuel_dt, dT_water_dt, dDO_dt]

    y0 = [T_ambient, O2_frac0, T_ambient, 0.0, T_ambient, DO_sat(T_ambient)]
    t_eval = np.linspace(0, t_max, 600)
    sol = solve_ivp(rhs, [0, t_max], y0, method='RK45', max_step=5.0, t_eval=t_eval)

    # 用 t_eval 重算 HRR(避免用 rhs 呼叫時的雜訊點)
    HRR_arr = []
    fuel_burned_arr = sol.y[3]
    O2_arr = sol.y[1]
    for i, t in enumerate(sol.t):
        HRR_growth = alpha_fire * t**2
        fuel_factor = 1.0 if fuel_burned_arr[i] < m_fuel_kg else 0.0
        O2_factor = np.clip((O2_arr[i] - extinct_lo) / (extinct_hi - extinct_lo), 0, 1)
        HRR_arr.append(min(HRR_growth * fuel_factor * O2_factor, 8000.0))

    return {
        "t": sol.t, "T_room": sol.y[0], "O2_frac": sol.y[1], "T_wall": sol.y[2],
        "fuel_burned": sol.y[3], "T_water": sol.y[4], "DO": sol.y[5], "HRR": np.array(HRR_arr)
    }

# 三個代表案例
cases = {
    "Clear survive\n(V=50 m3, R=3.0m)": simulate_full(50, 3.0),
    "Near boundary\n(V=150 m3, R=0.9m)": simulate_full(150, 0.9),
    "Clear death\n(V=250 m3, R=0.5m)": simulate_full(250, 0.5),
}

fig, axes = plt.subplots(5, 3, figsize=(16, 14), sharex='col')

for col, (title, r) in enumerate(cases.items()):
    t_min = r["t"] / 60
    axes[0, col].plot(t_min, r["HRR"], color='darkred')
    axes[0, col].set_title(title, fontsize=11)
    axes[0, col].set_ylabel("HRR (kW)")

    axes[1, col].plot(t_min, r["O2_frac"]*100, color='steelblue')
    axes[1, col].axhline(15, color='gray', ls='--', lw=1, label='extinction ~15%')
    axes[1, col].set_ylabel("Room O2 (%)")

    axes[2, col].plot(t_min, r["T_room"], color='orange', label='T_room')
    axes[2, col].plot(t_min, r["T_wall"], color='brown', ls=':', label='T_wall')
    axes[2, col].set_ylabel("Temp (C)")
    axes[2, col].legend(fontsize=8)

    axes[3, col].plot(t_min, r["T_water"], color='crimson')
    axes[3, col].axhline(T_lethal, color='red', ls='--', lw=1, label=f'lethal {T_lethal}C')
    axes[3, col].set_ylabel("T_water (C)")
    axes[3, col].legend(fontsize=8)

    axes[4, col].plot(t_min, r["DO"], color='green')
    axes[4, col].axhline(DO_lethal, color='red', ls='--', lw=1, label=f'lethal {DO_lethal}mg/L')
    axes[4, col].set_ylabel("DO (mg/L)")
    axes[4, col].set_xlabel("Time (min)")
    axes[4, col].legend(fontsize=8)

plt.suptitle("Time-history comparison: survive / boundary / death cases", fontsize=14)
plt.tight_layout()
plt.savefig("timehistory.png", dpi=140)
print("saved")

# 印出關鍵時間點供文字說明使用
for title, r in cases.items():
    max_Tw_idx = np.argmax(r["T_water"])
    min_DO_idx = np.argmin(r["DO"])
    print(f"\n=== {title.splitlines()[0]} ===")
    print(f"  max T_water = {r['T_water'].max():.1f} C at t={r['t'][max_Tw_idx]/60:.1f} min")
    print(f"  min DO      = {r['DO'].min():.2f} mg/L at t={r['t'][min_DO_idx]/60:.1f} min")
    ext_idx = np.argmax(r["O2_frac"] <= 0.15) if (r["O2_frac"]<=0.15).any() else None
    if ext_idx:
        print(f"  O2 hits 15% at t={r['t'][ext_idx]/60:.1f} min")
    else:
        print(f"  O2 never drops to 15% (min O2={r['O2_frac'].min()*100:.1f}%)")
