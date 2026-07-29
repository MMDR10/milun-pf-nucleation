"""
Final analysis: Pre-seismic P_f decline rate vs earthquake magnitude
"""
import numpy as np

# Our findings
depths = ['P170', 'P360', 'P424', 'P464']
slopes = [-1.53e-7, -9.55e-8, -3.60e-8, -6.33e-8]  # MPa/s
r2 = [0.836, 0.638, 0.321, 0.725]
total_drop = [-0.000258, -0.000160, -0.000060, -0.000106]

print("P_f Pre-seismic Decline Rate vs Depth")
print("=" * 60)
print(f"{'Depth':>8s} {'Slope(MPa/hr)':>16s} {'R²':>8s} {'28min_drop(kPa)':>16s} {'Rate(mm_w/s)':>16s}")
print("-" * 64)

for i, d in enumerate(depths):
    mpa_hr = slopes[i] * 3600 * 1000  # MPa/s → kPa/hr
    mm_ws = slopes[i] * 1e6 / (9.81 * 1000) * 1000  # MPa/s → mm water/s
    total_kpa = total_drop[i] * 1000
    print(f"{d:>8s} {mpa_hr:>16.4f} {r2[i]:>8.3f} {total_kpa:>16.4f} {mm_ws:>16.4f}")

# Gradient analysis
print("\nP_f Gradient Change Analysis")
print("=" * 60)
# Hydrostatic pre-seismic
print(f"P170-P360 gradient: ~1.0 (hydrostatic)")
print(f"P424 gradient anomaly: Z = -16.88 post-seismic")
print(f"")
print(f"Phase space divergence (P424×P464): θ +72.5% (p=0.0018)")
print(f"Phase space divergence (P170×P464): θ +54.3% (p=0.0186)")
print(f"")
print(f"Co-seismic step: ALL depths Z < |0.4| (no co-seismic change)")
print(f"Post-seismic: CONTINUED decline (faster than pre-seismic)")
print(f"")
print(f"Hydrostatic check:")
print(f"  P170-P360: ~1010 m/m (perfect)")
print(f"  P424-P464: ~847 m/m (sub-hydrostatic)")
