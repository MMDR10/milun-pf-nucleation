#!/usr/bin/env python3
"""
Ô-HAT 200Hz — 第二階段：精細前兆演化分析 + Coseismic 響應定量
"""
import numpy as np
import csv, json, os

DATA_FILE = '/tmp/pf_data.csv'
OUTPUT_DIR = '/app/working/workspaces/tygtDc/data/Milun_Pf'

with open(DATA_FILE, 'r') as f:
    reader = csv.reader(f)
    next(reader)
    rows = list(reader)

def parse_time(t_str):
    t_str = t_str.strip()
    if not t_str:
        return -1
    parts = t_str.split()
    hms = parts[0].split(':')
    h = int(hms[0]); m = int(hms[1]); s = int(hms[2]) if len(hms)>2 else 0
    if len(parts) > 1:
        if parts[1] == 'PM' and h < 12: h += 12
        elif parts[1] == 'AM' and h == 12: h = 0
    return h*3600 + m*60 + s

times = []; p424 = []; p464 = []; p170 = []; p360 = []
for r in rows:
    t = parse_time(r[0])
    if t < 0: continue
    times.append(t)
    p170.append(float(r[1])); p360.append(float(r[2]))
    p424.append(float(r[3])); p464.append(float(r[4]))

times = np.array(times); p424 = np.array(p424); p464 = np.array(p464)
p170 = np.array(p170); p360 = np.array(p360)

EQ_TIME = 7*3600 + 58*60 + 11  # 28691s
dt = 1/200  # 5ms

# =====================================================================
# ANALYSIS 1: Running θ spread with finer time resolution
# Use smaller windows: 30s window, 5s step
# =====================================================================
print("=== Time-resolved θ spread evolution ===")

window_s = 30
step_s = 5
win_samp = int(window_s * 200)
step_samp = int(step_s * 200)

def running_ohat(p, times, eq_t):
    """Compute running Ô-HAT metrics with fine time resolution"""
    n = len(p)
    dp = np.gradient(p, dt)
    d2p = np.gradient(dp, dt)
    
    results = []
    for i in range(0, n - win_samp, step_samp):
        tw = times[i:i+win_samp]
        dpw = dp[i:i+win_samp]
        d2pw = d2p[i:i+win_samp]
        
        theta = np.arctan2(d2pw, dpw)
        
        # θ spread: circular standard deviation
        theta_mean_angle = np.arctan2(np.mean(np.sin(theta)), np.mean(np.cos(theta)))
        theta_spread = np.sqrt(-2 * np.log(np.abs(np.mean(np.exp(1j * theta)))))
        
        # TopSV
        pw = p[i:i+win_samp]
        topsv = np.mean(np.abs(pw - np.mean(pw)))
        
        results.append({
            't_center': np.mean(tw),
            'theta_spread': float(theta_spread),
            'theta_mean': float(theta_mean_angle),
            'TopSV': float(topsv),
            'P_mean': float(np.mean(pw))
        })
    
    return results

p424_30s = running_ohat(p424, times, EQ_TIME)
p464_30s = running_ohat(p464, times, EQ_TIME)

# Convert to arrays for analysis
def to_array(results):
    return {k: np.array([r[k] for r in results]) for k in results[0].keys()}

arr424 = to_array(p424_30s)
arr464 = to_array(p464_30s)

# =====================================================================
# ANALYSIS 2: Pre-seismic θ spread anomaly — find onset time
# =====================================================================
print("\n=== Onset time detection ===")

# Baseline: 7:30 - 7:38 (first 8 min)
baseline_mask = arr424['t_center'] < EQ_TIME - 20*60
baseline_spread = np.mean(arr424['theta_spread'][baseline_mask])
baseline_std = np.std(arr424['theta_spread'][baseline_mask])
threshold = baseline_spread + 2*baseline_std  # 2-sigma threshold

print(f"P424 θ_spread baseline: {baseline_spread:.4f} ± {baseline_std:.4f}")
print(f"Threshold (2σ): {threshold:.4f}")

# Find when θ spread first exceeds threshold in the last 20 min
pre_mask = (arr424['t_center'] >= EQ_TIME - 20*60) & (arr424['t_center'] < EQ_TIME)
pre_times = arr424['t_center'][pre_mask]
pre_spread = arr424['theta_spread'][pre_mask]

exceed_thresh = pre_spread > threshold
if np.any(exceed_thresh):
    first_exceed_idx = np.where(exceed_thresh)[0][0]
    onset_time = pre_times[first_exceed_idx]
    onset_rel = EQ_TIME - onset_time
    print(f"P424 θ_spread 2σ exceed onset: T-{onset_rel:.0f}s (at {onset_time/3600:.1f}h)")
    print(f"  = {onset_rel/60:.1f} min before earthquake")
    
    # Max θ spread in pre-seismic
    max_idx = np.argmax(pre_spread)
    max_spread = pre_spread[max_idx]
    max_time_rel = EQ_TIME - pre_times[max_idx]
    print(f"Peak θ_spread: {max_spread:.4f} at T-{max_time_rel:.0f}s (z={(max_spread-baseline_spread)/baseline_std:.1f}σ)")
else:
    print("No 2σ exceedance found")

# Compare last 5 min vs early baseline
last5_mask = (arr424['t_center'] >= EQ_TIME - 5*60) & (arr424['t_center'] < EQ_TIME)
last5_spread = np.mean(arr424['theta_spread'][last5_mask])
last5_z = (last5_spread - baseline_spread) / baseline_std
print(f"\nP424 last 5 min mean spread: {last5_spread:.4f} (z={last5_z:.1f}σ)")

# Same for P464
baseline_mask464 = arr464['t_center'] < EQ_TIME - 20*60
base_spread464 = np.mean(arr464['theta_spread'][baseline_mask464])
base_std464 = np.std(arr464['theta_spread'][baseline_mask464])
last5_mask464 = (arr464['t_center'] >= EQ_TIME - 5*60) & (arr464['t_center'] < EQ_TIME)
last5_spread464 = np.mean(arr464['theta_spread'][last5_mask464])
last5_z464 = (last5_spread464 - base_spread464) / base_std464
print(f"P464 last 5 min mean spread: {last5_spread464:.4f} (z={last5_z464:.1f}σ)")

# =====================================================================
# ANALYSIS 3: Fine-grained θ spread ratio comparison  
# =====================================================================
print("\n=== Fine-grained θ spread ratio (pre-20min / baseline) ===")
for name, arr in [('P170', None), ('P360', None), ('P424', arr424), ('P464', arr464)]:
    if arr is None:
        continue
    bm = arr['theta_spread'][arr['t_center'] < EQ_TIME - 20*60]
    pm = arr['theta_spread'][(arr['t_center'] >= EQ_TIME - 20*60) & (arr['t_center'] < EQ_TIME)]
    if len(bm) > 0 and len(pm) > 0:
        ratio = np.mean(pm) / np.mean(bm)
        print(f"  {name}: baseline_spread={np.mean(bm):.4f}, pre_spread={np.mean(pm):.4f}, ratio={ratio:.3f}")

# =====================================================================
# ANALYSIS 4: Coseismic detailed window (07:58-08:02)
# =====================================================================
print("\n=== Coseismic Detail (07:58:00 - 08:02:00) ===")
t_start_cs = 7*3600 + 58*60  # 28680
t_end_cs = 7*3600 + 62*60    # 28920

cs_mask = (times >= t_start_cs) & (times <= t_end_cs)
cs_t = times[cs_mask] - EQ_TIME  # relative to EQ

for name, p in [('P170', p170), ('P360', p360), ('P424', p424), ('P464', p464)]:
    cs_p = p[cs_mask]
    pre_eq = np.mean(cs_p[cs_t < -5])  # 5s before EQ
    post_eq = np.mean(cs_p[cs_t > 10])  # 10s after EQ start
    min_p = np.min(cs_p)
    max_p = np.max(cs_p)
    drop = pre_eq - min_p
    print(f"  {name}: pre={pre_eq*100:.3f}bar, min={min_p*100:.3f}bar, max={max_p*100:.3f}bar, range={(max_p-min_p)*100:.2f}bar, drop={drop*100:.3f}bar")

# =====================================================================
# ANALYSIS 5: Pre-seismic P_f trend — linear regression on last 5 min
# =====================================================================
print("\n=== Pre-seismic P_f slope (last 5 min) ===")
ps_mask = (times >= EQ_TIME - 300) & (times < EQ_TIME)
ps_t = times[ps_mask] - EQ_TIME  # relative seconds
for name, p in [('P170', p170), ('P360', p360), ('P424', p424), ('P464', p464)]:
    ps_p = p[ps_mask]
    # Linear fit
    A = np.vstack([ps_t, np.ones_like(ps_t)]).T
    m, c = np.linalg.lstsq(A, ps_p, rcond=None)[0]
    # Residual
    residual = ps_p - (m*ps_t + c)
    rmse = np.sqrt(np.mean(residual**2))
    print(f"  {name}: slope={m*100:.6f} bar/s = {m*60*100:.4f} bar/min | RMSE={rmse*100:.4f}bar")

# =====================================================================
# ANALYSIS 6: Post-seismic vs pre-seismic P_f level change
# =====================================================================
print("\n=== Post-seismic P_f level change ===")
# Pre: last 60s before EQ
pre_60_mask = (times >= EQ_TIME - 60) & (times < EQ_TIME)
# Post: 60-120s after EQ  
post_60_mask = (times >= EQ_TIME + 60) & (times < EQ_TIME + 120)
for name, p in [('P170', p170), ('P360', p360), ('P424', p424), ('P464', p464)]:
    pre = np.mean(p[pre_60_mask])
    post = np.mean(p[post_60_mask])
    delta = post - pre
    print(f"  {name}: pre={pre*100:.4f}bar, post={post*100:.4f}bar, Δ={delta*100:.4f}bar ({(delta/pre)*100:.3f}%)")

print("\nDone!")
