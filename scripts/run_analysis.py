"""
Milun Fault P_f + Acceleration Ô-HAT Analysis
===============================================
Compute:
 1. P_f timeseries at 4 depths - basic stats + anomaly detection
 2. Cross-depth coupling (θ matrix)
 3. Fluid-Solid coupling: P_f vs Acceleration (at 392m depth since closest to P_f sensors)
 4. Pre-seismic anomaly at each depth
"""
import numpy as np
import pandas as pd
from scipy import signal, stats
import json

print("=" * 60)
print("MILUN FAULT 200Hz P_f ANALYSIS - 2024/04/03 M7.4 Hualien")
print("=" * 60)

# ── Load P_f data ──
base = '/app/working/workspaces/tygtDc/data/Milun_Pf/The 200 Hz Sampling Interval of Porewater Pressure'
df_p = pd.read_csv(f'{base}/Groundwater Pressure Data .csv')
df_a = pd.read_csv(f'{base}/Acceleration Data.csv')

# Rename columns for cleanliness
df_p.columns = ['Time', 'P170', 'P360', 'P424', 'P464']
df_a = df_a.rename(columns={
    'N-S Acceleration (𝑚∕𝑠^2 )(392 m)': 'Acc_NS_392',
    'E-W Acceleration (𝑚∕𝑠^2 ) (392 m)': 'Acc_EW_392',
    'Z Acceleration (𝑚∕𝑠^2 ) (392 m)': 'Acc_Z_392',
})

# ── Basic stats ──
print(f"\n📊 P_f Basic Statistics (MPa)")
print(f"{'Depth':>8s} {'Mean':>12s} {'Std':>12s} {'Min':>12s} {'Max':>12s} {'Range':>12s}")
print("-" * 68)
depths = ['P170', 'P360', 'P424', 'P464']
depth_m = [170, 360, 424, 464]
for d in depths:
    col = df_p[d]
    print(f"{d:>8s} {col.mean():>12.6f} {col.std():>12.6f} {col.min():>12.6f} {col.max():>12.6f} {(col.max()-col.min()):>12.6f}")

# ── Hydrostatic gradient check ──
print(f"\n📐 Hydrostatic Gradient Check")
for i in range(len(depths)-1):
    dP = df_p[depths[i+1]].mean() - df_p[depths[i]].mean()
    dz = depth_m[i+1] - depth_m[i]
    grad = dP / dz * 1e6 / 9.81  # Convert MPa/m to m water/m (rho*g)
    print(f"  {depths[i]}→{depths[i+1]}: dP={dP:.4f} MPa, dz={dz}m, grad={grad:.2f} m/m (hydrostatic=1.0)")

# ── Time array ──
# 720000 rows at 200 Hz = 3600 seconds from 7:30 AM
t_total = 3600.0  # seconds
t_arr = np.linspace(0, t_total, len(df_p))
eq_time = 28 * 60  # 7:58 AM = 28 min after start = 1680 seconds

# ── Anomaly detection: pre-seismic vs post-seismic ──
pre_mask = t_arr < eq_time
post_mask = t_arr >= eq_time

print(f"\n🔍 Pre-seismic (7:30-7:58) vs Co/Post-seismic (7:58-8:30) Anomalies")
print(f"{'Depth':>8s} {'Pre_mean':>12s} {'Post_mean':>12s} {'Diff':>12s} {'Pre_std':>12s} {'Z-score':>12s}")
print("-" * 68)
for d in depths:
    pre = df_p[d][pre_mask]
    post = df_p[d][post_mask]
    diff = post.mean() - pre.mean()
    z = diff / pre.std()
    print(f"{d:>8s} {pre.mean():>12.6f} {post.mean():>12.6f} {diff:>12.6f} {pre.std():>12.6f} {z:>12.4f}")

# ── Ô-HAT: intra-P_f cross-depth correlation ──
print(f"\n🔄 Ô-HAT Cross-Depth Correlation Matrix (r)")
print(f"{'':>8s}", end="")
for d in depths:
    print(f"{d:>12s}", end="")
print()
for d1 in depths:
    print(f"{d1:>8s}", end="")
    for d2 in depths:
        r, _ = stats.pearsonr(df_p[d1], df_p[d2])
        print(f"{r:>12.4f}", end="")
    print()

# ── Downsampled correlation (1 Hz decimated) for robustness ──
print(f"\n🔄 Ô-HAT Cross-Depth @ 1 Hz decimated (r)")
df_1hz = df_p.iloc[::200].reset_index(drop=True)
for d1 in depths:
    print(f"{d1:>8s}", end="")
    for d2 in depths:
        r, _ = stats.pearsonr(df_1hz[d1], df_1hz[d2])
        print(f"{r:>12.4f}", end="")
    print()

# ── Ô-HAT: Fluid-Solid coupling (P_f vs Acceleration around EQ) ──
print(f"\n⚡ FLUID-SOLID COUPLING: P_f vs Acceleration (at 392m)")
# Acceleration is 24000 rows = 120 sec = exactly at earthquake time
# Find the matching P_f window: 200 Hz * 120 sec = 24000 rows
# Acceleration starts at 07:58:00, P_f starts at 7:30:00
# So the P_f window is rows [eq_time*200 : eq_time*200 + 24000]
pf_start = int(eq_time * 200)
pf_end = pf_start + 24000

print(f"  P_f window: rows {pf_start}-{pf_end} ({pf_end-pf_start})")
print(f"  Acc window: rows 0-{len(df_a)} ({len(df_a)})")

p_f_window = df_p.iloc[pf_start:pf_end]

for depth in depths:
    print(f"\n  --- {depth} vs Acceleration ---")
    for acc_col in ['Acc_NS_392', 'Acc_EW_392', 'Acc_Z_392']:
        r, p_val = stats.pearsonr(p_f_window[depth][:len(df_a)], df_a[acc_col])
        print(f"    {depth} × {acc_col}: r={r:.4f}, p={p_val:.6f}")

# ── Pre-seismic sliding window anomaly ──
print(f"\n📉 PRE-SEISMIC SLIDING WINDOW (1-min windows)")
print(f"  Comparing each 1-min window against the full pre-seismic baseline")
print(f"{'Depth':>8s} {'Window':>14s} {'Mean':>12s} {'Baseline':>12s} {'Z':>10s} {'Anomaly?':>10s}")
print("-" * 68)

pre_data = df_p.iloc[:pf_start]
for depth in depths:
    baseline_mean = pre_data[depth].mean()
    baseline_std = pre_data[depth].std()
    
    # 1-min windows = 12000 rows
    n_windows = pf_start // 12000
    for w in range(n_windows):
        w_start = w * 12000
        w_end = w_start + 12000
        w_mean = pre_data[depth].iloc[w_start:w_end].mean()
        z = (w_mean - baseline_mean) / (baseline_std / np.sqrt(12000))
        minute = 30 + w  # 7:30 + w min
        anomaly = "🟢" if abs(z) < 3 else "🔴" if abs(z) < 5 else "🔥🔥"
        if abs(z) > 3:
            print(f"{depth:>8s} 7:{minute:02d}-{minute+1:02d}  {w_mean:>12.6f} {baseline_mean:>12.6f} {z:>10.2f} {anomaly:>10s}")

print(f"\n✅ Analysis complete!")
