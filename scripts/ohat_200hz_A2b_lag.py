#!/usr/bin/env python3
"""Fix: Run only the lagged cross-correlation (A2b)"""
import numpy as np
from scipy import stats, signal
import csv, json, os, time

DATA_FILE = '/tmp/pf_data.csv'
OUTPUT_DIR = '/app/working/workspaces/tygtDc/data/Milun_Pf'

with open(DATA_FILE, 'r') as f:
    reader = csv.reader(f)
    next(reader)
    rows = list(reader)

def parse_time(t_str):
    t_str = t_str.strip()
    if not t_str: return -1
    parts = t_str.split()
    hms = parts[0].split(':')
    h = int(hms[0]); m = int(hms[1]); s = int(hms[2]) if len(hms)>2 else 0
    if len(parts) > 1:
        if parts[1] == 'PM' and h < 12: h += 12
        elif parts[1] == 'AM' and h == 12: h = 0
    return h*3600 + m*60 + s

times = []; p424 = []; p464 = []
for r in rows:
    t = parse_time(r[0])
    if t < 0: continue
    times.append(t)
    p424.append(float(r[3])); p464.append(float(r[4]))

times = np.array(times); p424 = np.array(p424); p464 = np.array(p464)
EQ_TIME = 28691
fs = 200.0

pre_mask = times < EQ_TIME
post_mask = times >= EQ_TIME + 30

print("=== A2b: LAGGED CROSS-CORRELATION ===")

# Use decimated data for efficiency (200Hz → 2Hz)
dec = 100
p424_d = signal.decimate(p424, dec, ftype='fir')
p464_d = signal.decimate(p464, dec, ftype='fir')
times_d = times[::dec]
fs_d = 2.0

pre_d = pre_mask[::dec]
post_d = post_mask[::dec]

print(f"Decimated: {len(p424_d)} samples @ {fs_d} Hz")

# ===== Full pre-seismic lagged xcorr =====
max_lag_s = 120  # ±2 min
max_lag = int(max_lag_s * fs_d)

x_pre = p424_d[pre_d]; y_pre = p464_d[pre_d]
x_pre_n = (x_pre - np.mean(x_pre)) / np.std(x_pre)
y_pre_n = (y_pre - np.mean(y_pre)) / np.std(y_pre)

# Use only first 10 min of post for fair comparison
post_10min = int(10 * 60 * fs_d)
x_post = p424_d[post_d][:post_10min]; y_post = p464_d[post_d][:post_10min]
x_post_n = (x_post - np.mean(x_post)) / np.std(x_post)
y_post_n = (y_post - np.mean(y_post)) / np.std(y_post)

# Cross-correlation
cc_pre = np.correlate(x_pre_n, y_pre_n, mode='same') / len(x_pre_n)
cc_post = np.correlate(x_post_n, y_post_n, mode='same') / len(x_post_n)
lags = np.arange(-len(cc_pre)//2, len(cc_pre)//2 + 1)

# Extract lags of interest - ensure same length
n_lags = len(lags)
n_cc = len(cc_pre)
if n_lags > n_cc:
    lags = lags[:n_cc]
elif n_cc > n_lags:
    cc_pre = cc_pre[:n_lags]

lag_mask = np.abs(lags) <= max_lag
lags_zoom = lags[lag_mask]
cc_pre_zoom = cc_pre[lag_mask] / fs_d  # convert to seconds
cc_post_zoom_full = np.correlate(x_post_n, y_post_n, mode='same') / len(x_post_n)
lags_post = np.arange(-len(cc_post_zoom_full)//2, len(cc_post_zoom_full)//2 + 1)
n_lags_p = len(lags_post)
n_cc_p = len(cc_post_zoom_full)
if n_lags_p > n_cc_p:
    lags_post = lags_post[:n_cc_p]
elif n_cc_p > n_lags_p:
    cc_post_zoom_full = cc_post_zoom_full[:n_lags_p]
post_lag_mask = np.abs(lags_post) <= max_lag
cc_post_zoom = cc_post_zoom_full[post_lag_mask] / fs_d

print(f"\nPre-seismic lagged xcorr (full 28 min):")
print(f"  Peak r = {np.max(cc_pre_zoom):+.4f} @ lag = {lags_zoom[np.argmax(cc_pre_zoom)]:.0f}s")
print(f"  Min r  = {np.min(cc_pre_zoom):+.4f} @ lag = {lags_zoom[np.argmin(cc_pre_zoom)]:.0f}s")
r0_mask = lags_zoom == 0
if np.any(r0_mask):
    print(f"  r @ lag=0: {cc_pre_zoom[r0_mask][0]:+.4f}")
else:
    print(f"  r @ lag=0: N/A")

print(f"\nPost-seismic lagged xcorr (10 min):")
print(f"  Peak r = {np.max(cc_post_zoom):+.4f} @ lag = {lags_zoom[np.argmax(cc_post_zoom)]:.0f}s")
print(f"  Min r  = {np.min(cc_post_zoom):+.4f} @ lag = {lags_zoom[np.argmin(cc_post_zoom)]:.0f}s")
r0_mask_p = lags_zoom == 0
if np.any(r0_mask_p):
    print(f"  r @ lag=0: {cc_post_zoom[r0_mask_p][0]:+.4f}")
else:
    print(f"  r @ lag=0: N/A")

# ===== Also check last 5 min specifically =====
print(f"\nLast 5 min pre-EQ:")
last5 = (times >= EQ_TIME - 300) & (times < EQ_TIME)
x_l5 = p424[last5]; y_l5 = p464[last5]
x_l5_n = (x_l5 - np.mean(x_l5)) / np.std(x_l5)
y_l5_n = (y_l5 - np.mean(y_l5)) / np.std(y_l5)
max_lag_l5 = int(30 * fs)  # ±30s
cc_l5 = np.correlate(x_l5_n, y_l5_n, mode='same') / len(x_l5_n)
center_l5 = len(cc_l5)//2
lags_l5 = (np.arange(-center_l5, len(cc_l5)-center_l5)) / fs
zoom_l5 = np.abs(lags_l5) <= 30
print(f"  Peak r = {np.max(cc_l5[zoom_l5]):+.4f} @ lag = {lags_l5[zoom_l5][np.argmax(cc_l5[zoom_l5])]:.0f}s")
print(f"  r @ lag=0: {cc_l5[zoom_l5][lags_l5[zoom_l5]==0][0]:+.4f}")

# ===== Also check 2-min windows sliding lag =====
print(f"\n=== Sliding window lag analysis ===")
# For each 2-min window, find the lag that maximizes cross-correlation
win_samp = int(2*60*fs)
step_samp = int(30*fs)  # 30s step
pre_only = pre_mask

n_windows = (np.sum(pre_only) - win_samp) // step_samp + 1
print(f"  Analyzing {n_windows} pre-seismic windows...")

lag_evolution = []
for i in range(0, np.sum(pre_only) - win_samp, step_samp):
    xw = p424[pre_only][i:i+win_samp]
    yw = p464[pre_only][i:i+win_samp]
    tw = times[pre_only][i:i+win_samp]
    
    xw_n = (xw - np.mean(xw)) / np.std(xw)
    yw_n = (yw - np.mean(yw)) / np.std(yw)
    
    max_lag_w = int(30 * fs)  # ±30s
    cc_w = np.correlate(xw_n, yw_n, mode='same') / len(xw_n)
    center_w = len(cc_w)//2
    lags_w = (np.arange(-center_w, len(cc_w)-center_w)) / fs
    
    zoom = np.abs(lags_w) <= max_lag_w
    cc_z = cc_w[zoom]
    lags_z = lags_w[zoom]
    
    peak_idx = np.argmax(cc_z)
    lag_evolution.append({
        't_center': float(np.mean(tw)),
        'peak_r': float(cc_z[peak_idx]),
        'peak_lag_s': float(lags_z[peak_idx]),
        'r_at_lag0': float(cc_z[lags_z==0][0]) if np.any(lags_z==0) else 0
    })

# Compare early vs late pre-seismic
early_lag = [l for l in lag_evolution if l['t_center'] < EQ_TIME - 15*60]
late_lag = [l for l in lag_evolution if l['t_center'] >= EQ_TIME - 5*60]

if early_lag and late_lag:
    print(f"\n  Early pre-seismic (n={len(early_lag)}):")
    print(f"    Mean peak lag = {np.mean([l['peak_lag_s'] for l in early_lag]):.1f}s")
    print(f"    Mean r@lag0 = {np.mean([l['r_at_lag0'] for l in early_lag]):+.4f}")
    
    print(f"  Late pre-seismic (n={len(late_lag)}):")
    print(f"    Mean peak lag = {np.mean([l['peak_lag_s'] for l in late_lag]):.1f}s")
    print(f"    Mean r@lag0 = {np.mean([l['r_at_lag0'] for l in late_lag]):+.4f}")
    
    # Is there a phase shift toward the end?
    early_lag_vals = np.array([l['peak_lag_s'] for l in early_lag])
    late_lag_vals = np.array([l['peak_lag_s'] for l in late_lag])
    print(f"    Δpeak_lag (late-early) = {np.mean(late_lag_vals) - np.mean(early_lag_vals):+.1f}s")
    
    # Statistical test
    t_stat, p_lag = stats.ttest_ind(early_lag_vals, late_lag_vals)
    print(f"    t-test early vs late: p={p_lag:.4f} {'✅' if p_lag<0.05 else '❌'}")

print("\nDone!")
