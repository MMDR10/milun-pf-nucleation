#!/usr/bin/env python3
"""
Ô-HAT 200Hz P_f 高頻分析
目標：用 200 Hz P424/P464 孔隙水壓數據，進行高時間分辨率 Ô-HAT 相空間分析，
     捕捉 M7.4 (2024/04/03) 地震前 ~28 分鐘的前兆訊號結構。
"""
import numpy as np
import csv, time, os, json

DATA_FILE = '/tmp/pf_data.csv'
OUTPUT_DIR = '/app/working/workspaces/tygtDc/data/Milun_Pf'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =====================================================================
# STEP 1: 載入數據
# =====================================================================
print("=== Loading 200 Hz P_f data ===")
t0 = time.time()

with open(DATA_FILE, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

print(f"Loaded {len(rows)} rows in {time.time()-t0:.1f}s")
print(f"Columns: {header}")

# Parse time strings to seconds from midnight
def parse_time(t_str):
    t_str = t_str.strip()
    if not t_str or t_str == 'Time':
        return -1
    parts = t_str.split()
    hms = parts[0].split(':')
    h = int(hms[0])
    m = int(hms[1])
    s = int(hms[2]) if len(hms) > 2 else 0
    if len(parts) > 1:
        if parts[1] == 'PM' and h < 12:
            h += 12
        elif parts[1] == 'AM' and h == 12:
            h = 0
    return h * 3600 + m * 60 + s

# Parse data
times = []
p170 = []
p360 = []
p424 = []
p464 = []

for r in rows:
    t = parse_time(r[0])
    if t < 0:
        continue
    times.append(t)
    try:
        p170.append(float(r[1]))
        p360.append(float(r[2]))
        p424.append(float(r[3]))
        p464.append(float(r[4]))
    except (ValueError, IndexError):
        times.pop()
        continue

times = np.array(times)
p170 = np.array(p170)
p360 = np.array(p360)
p424 = np.array(p424)
p464 = np.array(p464)

print(f"\nParsed {len(times)} valid rows")
print(f"Time range: {times[0]}s to {times[-1]}s ({times[-1]-times[0]:.0f}s total)")
print(f"Sampling: {1/np.mean(np.diff(times)):.1f} Hz (target: 200 Hz)")

# M7.4 earthquake at 07:58:11 local time = 7*3600+58*60+11 = 28691s
EQ_TIME = 7*3600 + 58*60 + 11  # 28691s
eq_idx = np.argmin(np.abs(times - EQ_TIME))
print(f"Earthquake index: {eq_idx} (time={times[eq_idx]}s)")
print(f"Pre-seismic samples: {eq_idx}")
print(f"Pre-seismic duration: {times[eq_idx]-times[0]:.0f}s = {(times[eq_idx]-times[0])/60:.1f} min")

# =====================================================================
# STEP 2: Ô-HAT 相空間分析 — 滑動窗口
# =====================================================================
print("\n=== Ô-HAT Phase Space Analysis ===")

def compute_ohat(p_series, times_series, window_sec=60, step_sec=10):
    """
    Compute Ô-HAT metrics on sliding windows.
    For each window:
    1. Compute derivative dP/dt
    2. Build 2D phase space (P, dP/dt)
    3. Compute Ô = -divergence of phase space flow
    4. Compute θ = angle between consecutive phase space vectors
    """
    n = len(p_series)
    window_samples = int(window_sec * 200)  # 200 Hz
    step_samples = int(step_sec * 200)
    dt = np.mean(np.diff(times_series))
    
    # First derivative (central difference)
    dp = np.gradient(p_series, dt)
    # Second derivative
    d2p = np.gradient(dp, dt)
    
    print(f"Window: {window_sec}s ({window_samples} samples), Step: {step_sec}s ({step_samples} samples)")
    
    results = []
    n_windows = (n - window_samples) // step_samples + 1
    
    for i in range(0, n - window_samples, step_samples):
        # Window slices
        pw = p_series[i:i+window_samples]
        dpw = dp[i:i+window_samples]
        d2pw = d2p[i:i+window_samples]
        tw = times_series[i:i+window_samples]
        
        # Ô: magnitude of phase space vector divergence
        # Ô = sqrt(dP² + d²P²) — Laplacian magnitude in phase space
        O = np.mean(np.sqrt(dpw**2 + d2pw**2))
        
        # θ: direction of phase space flow
        # θ = arctan2(d²P, dP) → averaged
        theta = np.mean(np.arctan2(d2pw, dpw))
        
        # Entropy of phase space: spread of θ values
        theta_bins = np.linspace(-np.pi, np.pi, 37)  # 10° bins
        theta_hist = np.histogram(np.arctan2(d2pw, dpw), bins=theta_bins)[0]
        theta_hist = theta_hist / np.sum(theta_hist)
        theta_hist = theta_hist[theta_hist > 0]
        H_theta = -np.sum(theta_hist * np.log(theta_hist)) / np.log(36)  # Normalized entropy
        
        # TopSV: Topological Signature Vector
        # Mean absolute deviation of P around window mean
        TopSV = np.mean(np.abs(pw - np.mean(pw)))
        
        # Rate of change of TopSV (slope)
        if len(results) > 5:
            recent = results[-5:]
            TopSV_slope = (TopSV - np.mean([r['TopSV'] for r in recent])) / (step_sec * 5)
        else:
            TopSV_slope = 0
        
        # Mean P_f in window
        P_mean = np.mean(pw)
        
        results.append({
            't_center': np.mean(tw),
            't_start': tw[0],
            't_end': tw[-1],
            'O': float(O),
            'theta': float(theta),
            'H_theta': float(H_theta),
            'TopSV': float(TopSV),
            'TopSV_slope': float(TopSV_slope),
            'P_mean': float(P_mean),
            'P_std': float(np.std(pw))
        })
    
    return results

# Run Ô-HAT on all 4 depths
print("\nAnalyzing P424 (most sensitive depth from prior analysis)...")
p424_oh = compute_ohat(p424, times)

print("\nAnalyzing P464 (deeper compartment)...")
p464_oh = compute_ohat(p464, times)

print("\nAnalyzing P170...")
p170_oh = compute_ohat(p170, times)

print("\nAnalyzing P360...")
p360_oh = compute_ohat(p360, times)

# =====================================================================
# STEP 3: 前兆檢測
# =====================================================================
print("\n=== Pre-seismic Anomaly Detection ===")

def detect_anomaly(results, eq_time_sec, pre_window_min=20):
    """Detect pre-seismic anomaly by comparing recent vs baseline phase space"""
    pre_window_sec = pre_window_min * 60
    pre_start = eq_time_sec - pre_window_sec
    
    # Baseline: [-28 min, -20 min] before EQ
    baseline = [r for r in results if eq_time_sec - 28*60 <= r['t_center'] < eq_time_sec - pre_window_sec]
    # Pre-seismic window: [-20 min, 0]
    pres = [r for r in results if eq_time_sec - pre_window_sec <= r['t_center'] < eq_time_sec]
    
    if not baseline or not pres:
        return None
    
    keys = ['O', 'theta', 'H_theta', 'TopSV', 'P_mean', 'P_std']
    output = {}
    
    for k in keys:
        base_vals = np.array([r[k] for r in baseline])
        pre_vals = np.array([r[k] for r in pres])
        
        base_mean = np.mean(base_vals)
        base_std = np.std(base_vals)
        pre_mean = np.mean(pre_vals)
        
        if base_std > 0:
            z = (pre_mean - base_mean) / base_std
        else:
            z = 0
        
        pct_change = (pre_mean - base_mean) / abs(base_mean) * 100 if base_mean != 0 else 0
        
        output[k] = {
            'baseline_mean': float(base_mean),
            'baseline_std': float(base_std),
            'pre_mean': float(pre_mean),
            'z_score': float(z),
            'pct_change': float(pct_change)
        }
    
    # θ diversity: compare baseline θ histogram vs pre-seismic θ histogram
    base_theta = np.array([r['theta'] for r in baseline])
    pre_theta = np.array([r['theta'] for r in pres])
    
    theta_spread_base = np.std(base_theta)
    theta_spread_pre = np.std(pre_theta)
    
    output['theta_spread_change'] = {
        'baseline_std': float(theta_spread_base),
        'pre_std': float(theta_spread_pre),
        'ratio': float(theta_spread_pre / theta_spread_base if theta_spread_base > 0 else 0)
    }
    
    return output

# Run anomaly detection for all depths
for name, oh_data in [('P170', p170_oh), ('P360', p360_oh), ('P424', p424_oh), ('P464', p464_oh)]:
    print(f"\n=== {name} Pre-seismic Anomaly (window: 60s pre-EQ vs baseline) ===")
    anomaly = detect_anomaly(oh_data, EQ_TIME)
    if anomaly is None:
        print("  Insufficient data")
        continue
    for k, v in anomaly.items():
        if isinstance(v, dict):
            if 'z_score' in v:
                print(f"  {k:12s}: baseline={v['baseline_mean']:.3e} → pre={v['pre_mean']:.3e} | z={v['z_score']:.2f} | {v['pct_change']:.1f}%")
            else:
                print(f"  {k:20s}: baseline_std={v['baseline_std']:.4f} → pre_std={v['pre_std']:.4f} | ratio={v['ratio']:.2f}")

# =====================================================================
# STEP 4: High-resolution time evolution plot (saved as JSON)
# =====================================================================
print("\n=== Saving results ===")

output = {
    'meta': {
        'data': '200 Hz P_f from Milun Fault',
        'source': 'Mendeley dn4338tv9g Marginingsih et al. 2026',
        'analysis': 'Ô-HAT sliding window (60s window, 10s step)',
        'eq_time': EQ_TIME,
        'eq_time_str': '2024-04-03 07:58:11 local'
    },
    'results': {}
}

for name, oh_data in [('P170', p170_oh), ('P360', p360_oh), ('P424', p424_oh), ('P464', p464_oh)]:
    anomaly = detect_anomaly(oh_data, EQ_TIME)
    output['results'][name] = {
        'timeseries': oh_data,
        'anomaly': anomaly
    }

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return super().default(obj)

with open(os.path.join(OUTPUT_DIR, 'ohat_200hz_results.json'), 'w') as f:
    json.dump(output, f, indent=2, cls=NumpyEncoder)

print(f"Results saved to {OUTPUT_DIR}/ohat_200hz_results.json")

# =====================================================================
# STEP 5: Key findings summary
# =====================================================================
print("\n" + "="*60)
print("KEY FINDINGS")
print("="*60)

# 1. P424 theta spread ratio
for name in ['P424', 'P464']:
    a = output['results'][name]['anomaly']
    if a:
        tsr = a.get('theta_spread_change', {}).get('ratio', 0)
        top_z = a.get('TopSV', {}).get('z_score', 0)
        print(f"\n{name}:")
        print(f"  θ spread ratio (pre/baseline): {tsr:.2f}")
        print(f"  TopSV z-score: {top_z:.2f}")
        print(f"  P_mean change: {a['P_mean']['pct_change']:.2f}% (z={a['P_mean']['z_score']:.2f})")
        print(f"  Ô change: {a['O']['pct_change']:.2f}% (z={a['O']['z_score']:.2f})")
        print(f"  θ entropy change: {a['H_theta']['pct_change']:.2f}% (z={a['H_theta']['z_score']:.2f})")

print("\nDone!")
