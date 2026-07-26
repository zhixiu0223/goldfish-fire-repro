import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fish_ode_model2 import find_R_crit

V_room_range = np.linspace(20, 300, 12)

baseline = dict(chi_r=0.30, h_conv=20.0, kLa_mult=1.0)

# 每個係數測 0.5x / 1x(baseline) / 2x
param_tests = {
    "chi_r (輻射熱分率)": [
        ("0.15 (半)", dict(chi_r=0.15, h_conv=20.0, kLa_mult=1.0)),
        ("0.30 (基準)", baseline),
        ("0.60 (倍)", dict(chi_r=0.60, h_conv=20.0, kLa_mult=1.0)),
    ],
    "h_conv (對流係數)": [
        ("10 (半)", dict(chi_r=0.30, h_conv=10.0, kLa_mult=1.0)),
        ("20 (基準)", baseline),
        ("40 (倍)", dict(chi_r=0.30, h_conv=40.0, kLa_mult=1.0)),
    ],
    "kLa (氣液傳質)": [
        ("0.5x (慢)", dict(chi_r=0.30, h_conv=20.0, kLa_mult=0.5)),
        ("1.0x (基準)", baseline),
        ("2.0x (快)", dict(chi_r=0.30, h_conv=20.0, kLa_mult=2.0)),
    ],
}

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), sharey=True)

results_summary = {}
for ax, (pname, cases) in zip(axes, param_tests.items()):
    for label, params in cases:
        R_crits = [find_R_crit(V, params) for V in V_room_range]
        ax.plot(V_room_range, R_crits, marker='o', label=label)
        results_summary[(pname, label)] = R_crits
        print(f"{pname:20s} {label:12s} R_crit range: {min(R_crits):.2f} ~ {max(R_crits):.2f} m")
    ax.set_xlabel('Room volume (m^3)')
    ax.set_title(pname.split(' ')[0])
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

axes[0].set_ylabel('Critical distance R_crit (m)\n(above line = survive, below = death)')

plt.suptitle('Sensitivity of the death/survival boundary to chi_r, h_conv, kLa', fontsize=13)
plt.tight_layout()
plt.savefig('sensitivity.png', dpi=150)
print("\nplot saved")
