#!/usr/bin/env python3
"""
方向 A2 深度追 — P424×P464 交叉相關時變精細分析

目標：
1. 不同窗口大小對比（2min / 5min / 10min）
2. 統計顯著性檢驗（trend test on last 10 min r values）
3. 相位領先/落後分析（lagged cross-correlation）
4. 三階段劃分：early / mid / late pre-seismic
5. Change-point detection on r(t)
"""
import numpy as np
from scipy import stats, signal
import csv, json, os, time

DATA_FILE = '/tmp/pf_data.csv'
OUTPUT_DIR = '/app/working/workspaces/tygtDc/data/Milun_Pf'
os.makedirs(OUTPUT_DIR, exist_ok=True)

t0 = time.time()
print("=== Loading data ===")
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
EQ_TIME = 28691  # 07:58:11
fs = 200.0

# Masks
pre_mask = times < EQ_TIME
post_mask = times >= EQ_TIME + 30  # skip 30s co-seismic

# =====================================================================
# ANALYSIS 1: Multi-window rolling cross-correlation
# =====================================================================
print("\n" + "="*60)
print("ANALYSIS 1: Multi-window rolling cross-correlation")
print("="*60)

def rolling_r(x, y, times_arr, window_min, step_min=1, fs=200):
    """Rolling Pearson r with different window sizes"""
    win_samp = int(window_min * 60 * fs)
    step_samp = int(step_min * 60 * fs)
    n = min(len(x), len(y))
    
    results = []
    for i in range(0, n - win_samp, step_samp):
        xw = x[i:i+win_samp]
        yw = y[i:i+win_samp]
        r, p = stats.pearsonr(xw, yw)
        results.append({
            't_center': float(np.mean(times_arr[i:i+win_samp])),
            'r': float(r),
            'p': float(p),
            'z_r': float(np.arctanh(r))  # Fisher z-transform for normality
        })
    return results

for w_min in [2, 5, 10]:
    rr = rolling_r(p424, p464, times, window_min=w_min)
    
    # Pre vs post comparison
    pre_rr = [r for r in rr if r['t_center'] < EQ_TIME]
    post_rr = [r for r in rr if r['t_center'] > EQ_TIME + 60]
    
    if pre_rr and post_rr:
        pre_r_vals = np.array([r['r'] for r in pre_rr])
        post_r_vals = np.array([r['r'] for r in post_rr])
        
        # Pre-seismic three stages
        early = np.array([r['r'] for r in pre_rr if r['t_center'] < EQ_TIME - 18*60])
        mid = np.array([r['r'] for r in pre_rr if (r['t_center'] >= EQ_TIME - 18*60) & (r['t_center'] < EQ_TIME - 6*60)])
        late = np.array([r['r'] for r in pre_rr if r['t_center'] >= EQ_TIME - 6*60])
        
        print(f"\n  Window={w_min}min (n_pre={len(pre_r_vals)}, n_post={len(post_r_vals)}):")
        print(f"    Pre-seismic:       r̄={np.mean(pre_r_vals):+.4f} ± {np.std(pre_r_vals):.4f}")
        print(f"      Early (T-28~T-18): r̄={np.mean(early):+.4f} ± {np.std(early):.4f}" if len(early) > 0 else "")
        print(f"      Mid   (T-18~T-6):  r̄={np.mean(mid):+.4f} ± {np.std(mid):.4f}" if len(mid) > 0 else "")
        print(f"      Late  (T-6~T-0):   r̄={np.mean(late):+.4f} ± {np.std(late):.4f}" if len(late) > 0 else "")
        print(f"    Post-seismic:      r̄={np.mean(post_r_vals):+.4f} ± {np.std(post_r_vals):.4f}")
        print(f"    Δr (post-pre):     {np.mean(post_r_vals)-np.mean(pre_r_vals):+.4f}")
        
        # Mann-Whitney U test (pre vs post)
        if len(pre_r_vals) > 3 and len(post_r_vals) > 3:
            u_stat, p_val = stats.mannwhitneyu(pre_r_vals, post_r_vals, alternative='two-sided')
            print(f"    MWU test: U={u_stat:.0f}, p={p_val:.4f} {'✅' if p_val < 0.05 else '❌'}")
        
        # Trend test on late pre-seismic r values
        if len(late) > 3:
            late_times = np.array([r['t_center'] for r in pre_rr if r['t_center'] >= EQ_TIME - 6*60])
            late_r = np.array([r['r'] for r in pre_rr if r['t_center'] >= EQ_TIME - 6*60])
            slope, intercept, r_val, p_val_trend, std_err = stats.linregress(late_times, late_r)
            print(f"    Late trend: slope={slope*60:+.4f}/min, p={p_val_trend:.4f} {'✅' if p_val_trend < 0.05 else '❌'}")

# =====================================================================
# ANALYSIS 2: Lagged cross-correlation (pre vs post)
# =====================================================================
print("\n" + "="*60)
print("ANALYSIS 2: Lagged cross-correlation (P424 leads/lags P464)")
print("="*60)

def lagged_xcorr(x, y, max_lag_sec=60, fs=200):
    """Compute cross-correlation as function of lag"""
    max_lag = int(max_lag_sec * fs)
    # Use only pre-seismic data
    x_pre = x[pre_mask]
    y_pre = y[pre_mask]
    
    x_post = x[post_mask]
    y_post = y[post_mask]
    
    # Normalize
    x_pre_n = (x_pre - np.mean(x_pre)) / np.std(x_pre)
    y_pre_n = (y_pre - np.mean(y_pre)) / np.std(y_pre)
    x_post_n = (x_post - np.mean(x_post)) / np.std(x_post)
    y_post_n = (y_post - np.mean(y_post)) / np.std(y_post)
    
    # Cross-correlation via correlate
    lagged_r_pre = np.correlate(x_pre_n, y_pre_n, mode='same') / len(x_pre_n)
    lagged_r_post = np.correlate(x_post_n, y_post_n, mode='same') / len(x_post_n)
    lags = np.arange(-max_lag, max_lag + 1)
    
    center = len(lagged_r_pre) // 2
    lag_indices = np.arange(center - max_lag, center + max_lag + 1)
    
    print(f"  Pre-seismic peak lagged r: max={np.max(lagged_r_pre[lag_indices]):.4f} @ lag={lags[np.argmax(lagged_r_pre[lag_indices])]}s")
    print(f"  Post-seismic peak lagged r: max={np.max(lagged_r_post[lag_indices]):.4f} @ lag={lags[np.argmax(lagged_r_post[lag_indices])]}s")
    print(f"  Pre-seismic min lagged r: min={np.min(lagged_r_pre[lag_indices]):.4f} @ lag={lags[np.argmin(lagged_r_pre[lag_indices])]}s")
    
    # Also check last 5 min pre-seismic specifically for any phase shift
    last5_mask = (times >= EQ_TIME - 300) & (times < EQ_TIME)
    if np.sum(last5_mask) > fs * 10:  # at least 10s
        x_l5 = x[last5_mask]; y_l5 = y[last5_mask]
        x_l5_n = (x_l5 - np.mean(x_l5)) / np.std(x_l5)
        y_l5_n = (y_l5 - np.mean(y_l5)) / np.std(y_l5)
        max_lag_l5 = min(30, len(x_l5_n)//4)
        cc_l5 = np.correlate(x_l5_n, y_l5_n, mode='same') / len(x_l5_n)
        c_l5 = len(cc_l5)//2
        lags_l5 = np.arange(-max_lag_l5, max_lag_l5 + 1)
        cc_seg = cc_l5[c_l5 - max_lag_l5 : c_l5 + max_lag_l5 + 1]
        print(f"  Last 5 min pre-EQ peak lagged r: {np.max(cc_seg):.4f} @ lag={lags_l5[np.argmax(cc_seg)]}s")

# =====================================================================
# ANALYSIS 3: Coherence spectrum (pre vs post)
# =====================================================================
print("\n" + "="*60)
print("ANALYSIS 3: Magnitude-squared coherence")
print("="*60)

# Compute coherence between P424 and P464
f_pre, C_pre = signal.coherence(p424[pre_mask], p464[pre_mask], fs=fs, nperseg=fs*60)  # 60s windows
# Decimate post for speed (use first 10 min)
post_10min = int(10*60*fs)
f_post, C_post = signal.coherence(p424[post_mask][:post_10min], p464[post_mask][:post_10min], fs=fs, nperseg=fs*60)

print(f"  Mean coherence (0-10 Hz): pre={np.mean(C_pre[f_pre<=10]):.4f}, post={np.mean(C_post[f_post<=10]):.4f}")
print(f"  Mean coherence (10-50 Hz): pre={np.mean(C_pre[(f_pre>10)&(f_pre<=50)]):.4f}, post={np.mean(C_post[(f_post>10)&(f_post<=50)]):.4f}")

# Coherence peak comparison
print(f"  Max coherence pre: {np.max(C_pre):.4f} @ {f_pre[np.argmax(C_pre)]:.2f} Hz")
print(f"  Max coherence post: {np.max(C_post):.4f} @ {f_post[np.argmax(C_post)]:.2f} Hz")

# =====================================================================
# ANALYSIS 4: Statistical significance of the pre→post change
# =====================================================================
print("\n" + "="*60)
print("ANALYSIS 4: Bootstrap test for r change significance")
print("="*60)

def bootstrap_r_diff(x, y, pre_mask, post_mask, n_bootstrap=1000):
    """Bootstrap test: is the pre→post r change significant?"""
    x_pre, y_pre = x[pre_mask], y[pre_mask]
    x_post, y_post = x[post_mask], y[post_mask]
    
    # Observed difference
    r_pre, _ = stats.pearsonr(x_pre, y_pre)
    r_post, _ = stats.pearsonr(x_post, y_post[:len(x_post)])  # ensure same length if needed
    # For large data, use subsample for consistency
    n_sample = min(len(x_pre), len(x_post))
    np.random.seed(42)
    
    # Bootstrap: shuffle time labels
    diffs = []
    for _ in range(n_bootstrap):
        # Randomly assign each time point to pre or post
        idx = np.random.permutation(len(x_pre) + len(x_post))
        # Use shorter length for bootstrap efficiency
        n_short = min(len(x_pre), len(x_post), 10000)
        shuffled = np.concatenate([x_pre[:n_short], x_post[:n_short]])
        shuffled_y = np.concatenate([y_pre[:n_short], y_post[:n_short]])
        
        r1, _ = stats.pearsonr(shuffled[:n_short], shuffled_y[:n_short])
        r2, _ = stats.pearsonr(shuffled[n_short:], shuffled_y[n_short:])
        diffs.append(r2 - r1)
    
    diffs = np.array(diffs)
    obs_diff = r_post - r_pre
    p_boot = np.mean(np.abs(diffs) >= np.abs(obs_diff))
    
    return obs_diff, p_boot, np.percentile(diffs, [2.5, 97.5])

# Use decimated data for bootstrap (too large otherwise)
dec = 100  # decimate to 2 Hz
p424_d = signal.decimate(p424, dec, ftype='fir')
p464_d = signal.decimate(p464, dec, ftype='fir')
pre_mask_d = pre_mask[::dec]
post_mask_d = post_mask[::dec]

obs_diff, p_boot, ci = bootstrap_r_diff(p424_d, p464_d, pre_mask_d, post_mask_d, n_bootstrap=500)
print(f"  Observed Δr (post-pre): {obs_diff:+.4f}")
print(f"  Bootstrap 95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]")
print(f"  Bootstrap p-value: {p_boot:.4f} {'✅' if p_boot < 0.05 else '❌'}")

# =====================================================================
# ANALYSIS 5: Change-point detection in r(t) time series
# =====================================================================
print("\n" + "="*60)
print("ANALYSIS 5: Change-point detection in r(t)")
print("="*60)

# Use 2-min window for best time resolution
rr_2min = rolling_r(p424, p464, times, window_min=2, step_min=0.5)
rr_t = np.array([r['t_center'] for r in rr_2min])
rr_r = np.array([r['r'] for r in rr_2min])

# Simple change point: find when r crosses from near-zero to positive after EQ
# Look for the first time after EQ where r exceeds 0.3 (indicating onset of coupling)
post_rr = np.array([r for r in rr_2min if r['t_center'] > EQ_TIME])
if len(post_rr) > 0:
    post_t = np.array([r['t_center'] for r in post_rr])
    post_r = np.array([r['r'] for r in post_rr])
    
    # Find first time r > 0.3 with 3 consecutive points
    above_thresh = post_r > 0.3
    for i in range(len(above_thresh) - 3):
        if np.all(above_thresh[i:i+3]):
            coupling_onset = post_t[i]
            print(f"  Post-seismic coupling onset (r>0.3×3): T+{coupling_onset - EQ_TIME:.0f}s ({coupling_onset - EQ_TIME:.1f} min after EQ)")
            break
    else:
        print(f"  No sustained r>0.3 in post-seismic (max r={np.max(post_r):.3f})")

# Pre-seismic: find when r trends most negative
pre_rr_arr = np.array([r for r in rr_2min if r['t_center'] < EQ_TIME])
if len(pre_rr_arr) > 10:
    pre_t = np.array([r['t_center'] for r in pre_rr_arr])
    pre_r = np.array([r['r'] for r in pre_rr_arr])
    
    # Sliding window trend test
    for window_min in [3, 5]:
        n_win = int(window_min * 2)  # 2 points per minute (30s step)
        if n_win >= len(pre_r): continue
        
        best_z = 0
        best_t = 0
        for i in range(len(pre_r) - n_win):
            seg_t = pre_t[i:i+n_win]
            seg_r = pre_r[i:i+n_win]
            slope, _, _, p_val, _ = stats.linregress(seg_t, seg_r)
            # z-score: how many std devs from zero slope?
            if p_val < 0.1:  # marginally significant
                if abs(slope) > abs(best_z):
                    best_z = slope
                    best_t = np.mean(seg_t)
        
        if best_z != 0:
            direction = "declining" if best_z < 0 else "rising"
            print(f"  Most significant r-trend ({window_min}min window): {direction} at T-{EQ_TIME-best_t:.0f}s (z~{best_z*3600:.2f}/hr, p<0.10)")

# =====================================================================
# SAVE
# =====================================================================
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("\nDone!")
