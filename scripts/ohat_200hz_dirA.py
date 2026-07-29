#!/usr/bin/env python3
"""
Ô-HAT 200Hz — 方向 A 完整分析
A1: 帶通濾波分頻段
A2: P424×P464 交叉相關時變
A3: 固體潮/低頻 tidal response
A4: Spectrogram 高頻結構
"""
import numpy as np
from scipy import signal, stats
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
    p424.append(float(r[3]))
    p464.append(float(r[4]))

times = np.array(times); p424 = np.array(p424); p464 = np.array(p464)
EQ_TIME = 7*3600 + 58*60 + 11  # 28691s
dt = 1.0 / 200.0
fs = 200.0

# Pre-seismic mask
eq_idx = np.argmin(np.abs(times - EQ_TIME))
pre_mask = times < EQ_TIME
post_mask = times >= EQ_TIME

print(f"Loaded {len(times)} samples, {time.time()-t0:.1f}s")
print(f"Pre-seismic: {np.sum(pre_mask)} samples ({np.sum(pre_mask)/fs/60:.1f} min)")
print(f"Post-seismic: {np.sum(post_mask)} samples ({np.sum(post_mask)/fs/60:.1f} min)")

# =====================================================================
# A1: BANDPASS FILTER — 分頻段 Ô-HAT
# =====================================================================
print("\n" + "="*60)
print("A1: BANDPASS FILTER — Frequency-dependent θ spread")
print("="*60)

def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a

def bandpass_filter(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    return signal.filtfilt(b, a, data)

def band_ohat(p, times, fs, band_name, lowcut, highcut, win_sec=60, step_sec=10):
    """Compute θ spread ratio per frequency band
    
    Different frequency bands need different window sizes:
    - High freq (>1 Hz): 60s windows work fine
    - VLF (0.1-1 Hz): need 300s windows to capture low-frequency cycles
    - ULF (<0.1 Hz): need 600s+ windows (but we only have 28 min pre-seismic)
    """
    # Auto-select window based on lowest frequency
    if lowcut > 0:
        min_period = 1.0 / lowcut  # seconds
        # Window should be at least 5x the longest period
        recommended_win = max(win_sec, min_period * 5)
    else:
        recommended_win = win_sec
    
    # Cap at available data
    available = len(p) / fs
    recommended_win = min(recommended_win, available * 0.8)
    
    if recommended_win > win_sec * 1.5:
        print(f"    [{band_name}] Adjusting window: {win_sec}s → {recommended_win:.0f}s")
        win_sec_actual = recommended_win
    else:
        win_sec_actual = win_sec
    
    win_samp = int(win_sec_actual * fs)
    step_samp = int(step_sec * fs)
    n = len(p)
    
    if win_samp >= n:
        print(f"    [{band_name}] WARNING: window {win_samp} > data {n}, skipping")
        return []
    
    # Derivative
    dp = np.gradient(p, 1/fs)
    d2p = np.gradient(dp, 1/fs)
    
    results = []
    for i in range(0, n - win_samp, step_samp):
        dpw = dp[i:i+win_samp]
        d2pw = d2p[i:i+win_samp]
        theta = np.arctan2(d2pw, dpw)
        
        # Circular std of θ - handle NaN
        valid = ~np.isnan(theta)
        if np.sum(valid) < 10:
            continue
        theta_v = theta[valid]
        
        r_mean = np.mean(np.exp(1j * theta_v))
        if abs(r_mean) > 0:
            theta_spread = np.sqrt(-2 * np.log(abs(r_mean)))
        else:
            theta_spread = np.pi  # max spread
        
        # Mean Ô
        O = np.mean(np.sqrt(dpw[valid]**2 + d2pw[valid]**2))
        
        results.append({
            't_center': np.mean(times[i:i+win_samp]),
            'theta_spread': float(theta_spread),
            'O': float(O)
        })
    return results

# Bands to test
bands = [
    ('Wideband', 0.01, 95.0),     # Full range (DC and near-Nyquist removed)
    ('Tremor', 0.5, 5.0),         # Typical tremor band (volcanic/ tectonic tremor)
    ('LF', 1.0, 10.0),            # Low freq
    ('MF1', 10.0, 30.0),          # Mid freq 1
    ('MF2', 30.0, 60.0),          # Mid freq 2
    ('HF', 60.0, 95.0),           # High freq (near-Nyquist)
]

results_a1 = {}

for band_name, low, high in bands:
    # Only analyze pre-seismic (faster)
    p_pre = p424[pre_mask]
    t_pre = times[pre_mask]
    
    # Skip bands where filter can't be applied properly
    nyq = fs / 2  # 100 Hz
    if high >= nyq:
        if low > 0:
            # High cut at or above Nyquist: highpass only
            b, a = signal.butter(4, low / nyq, btype='high')
            p_filt = signal.filtfilt(b, a, p_pre)
        else:
            p_filt = p_pre.copy()
    elif low < 1e-5:
        print(f"  [{band_name}] Too low for filter, skipping")
        continue
    else:
        # Filter
        try:
            p_filt = bandpass_filter(p_pre, low, high, fs)
        except Exception as e:
            print(f"  [{band_name}] Filter failed: {e}, skipping")
            continue
    
    # Ô-HAT on filtered signal
    oh = band_ohat(p_filt, t_pre, fs, band_name, low, high)
    
    # Compute baseline vs pre-seismic ratio
    baseline = [r for r in oh if r['t_center'] < EQ_TIME - 20*60]
    pres = [r for r in oh if r['t_center'] >= EQ_TIME - 20*60]
    
    if baseline and pres:
        base_spread = np.mean([r['theta_spread'] for r in baseline])
        pre_spread = np.mean([r['theta_spread'] for r in pres])
        ratio = pre_spread / base_spread if base_spread > 0 else 0
        
        # Also check for enhanced spread in last 5 min
        last5 = [r for r in oh if r['t_center'] >= EQ_TIME - 5*60]
        last5_spread = np.mean([r['theta_spread'] for r in last5]) if last5 else 0
        last5_z = (last5_spread - base_spread) / (np.std([r['theta_spread'] for r in baseline]) if len(baseline) > 1 else 0.001)
        
        results_a1[band_name] = {
            'band_hz': f"{low}-{high}",
            'baseline_θ_spread': round(base_spread, 4),
            'pre-20min_θ_spread': round(pre_spread, 4),
            'θ_spread_ratio': round(ratio, 3),
            'last5min_θ_spread': round(last5_spread, 4),
            'last5min_z': round(last5_z, 2),
            'n_baseline': len(baseline),
            'n_pre': len(pres),
        }
        
        status = "✅" if ratio > 1.05 else ("⚠️" if ratio > 0.98 else "❌")
        print(f"  {band_name:10s} [{low:6.4f}-{high:5.1f} Hz]: "
              f"base={base_spread:.4f} → pre={pre_spread:.4f} | "
              f"ratio={ratio:.3f} {status} | last5min z={last5_z:+.2f}")

# =====================================================================
# A2: P424 × P464 ROLLING CROSS-CORRELATION
# =====================================================================
print("\n" + "="*60)
print("A2: P424 × P464 ROLLING CROSS-CORRELATION")
print("="*60)

def rolling_xcorr(x, y, window_min=5, step_min=1, fs=200):
    """Rolling Pearson correlation between two time series"""
    win_samp = int(window_min * 60 * fs)
    step_samp = int(step_min * 60 * fs)
    n = min(len(x), len(y))
    
    results = []
    for i in range(0, n - win_samp, step_samp):
        xw = x[i:i+win_samp]
        yw = y[i:i+win_samp]
        
        r, p = stats.pearsonr(xw, yw)
        
        results.append({
            't_center': np.mean(times[i:i+win_samp]),
            'r': float(r),
            'p': float(p),
            'r_abs': float(abs(r))
        })
    return results

# Full time range
xcorr = rolling_xcorr(p424, p464, window_min=5, step_min=1)

# Compare pre vs post for correlation change
pre_xc = [r for r in xcorr if r['t_center'] < EQ_TIME]
post_xc = [r for r in xcorr if r['t_center'] >= EQ_TIME + 60]  # skip 1 min co-seismic

if pre_xc and post_xc:
    pre_r = np.mean([r['r'] for r in pre_xc])
    post_r = np.mean([r['r'] for r in post_xc])
    pre_std = np.std([r['r'] for r in pre_xc])
    post_std = np.std([r['r'] for r in post_xc])
    
    print(f"  Pre-seismic mean r = {pre_r:.4f} ± {pre_std:.4f}")
    print(f"  Post-seismic mean r = {post_r:.4f} ± {post_std:.4f}")
    print(f"  Δr = {post_r - pre_r:+.4f}")
    
    # Check if correlation drops in last 10 min before EQ
    last10 = [r for r in pre_xc if r['t_center'] >= EQ_TIME - 10*60]
    if last10:
        last10_r = np.mean([r['r'] for r in last10])
        z = (last10_r - pre_r) / pre_std if pre_std > 0 else 0
        print(f"  Last 10 min pre-EQ mean r = {last10_r:.4f} (z = {z:.2f})")
    
    # Check if correlation increases after EQ (both sensors respond to shaking)
    early_post = [r for r in post_xc if r['t_center'] < EQ_TIME + 10*60]
    if early_post:
        early_post_r = np.mean([r['r'] for r in early_post])
        print(f"  Early post-seismic (10 min) mean r = {early_post_r:.4f}")

# =====================================================================
# A3: TIDAL/LOW-FREQUENCY RESPONSE
# =====================================================================
print("\n" + "="*60)
print("A3: LOW-FREQUENCY SPECTRAL ANALYSIS")
print("="*60)

# Use pre-seismic data only, decimate to 1 Hz for LF analysis
from scipy import signal as sp_signal

# Decimate to 1 Hz for computational efficiency
dec_factor = 200
p424_1hz = sp_signal.decimate(p424[pre_mask], dec_factor, ftype='fir')
p464_1hz = sp_signal.decimate(p464[pre_mask], dec_factor, ftype='fir')
n_1hz = len(p424_1hz)

# Welch PSD on 1 Hz data
f, psd_424 = sp_signal.welch(p424_1hz, fs=1, nperseg=min(600, n_1hz//4), noverlap=300)
f, psd_464 = sp_signal.welch(p464_1hz, fs=1, nperseg=min(600, n_1hz//4), noverlap=300)

print(f"  PSD computed: {len(f)} freq bins, {f[0]:.5f} - {f[-1]:.5f} Hz")

# Key tidal frequencies
tidal_freqs = {
    'M2 (main lunar)': 1.0 / (12.42 * 3600),    # 2.236e-5 Hz
    'S2 (main solar)': 1.0 / (12.0 * 3600),      # 2.315e-5 Hz
    'O1 (lunar diurnal)': 1.0 / (25.82 * 3600),  # 1.076e-5 Hz
    'K1 (lunisolar)': 1.0 / (23.93 * 3600),       # 1.161e-5 Hz
}

print("\n  Tidal frequency check (need >~1 hr data to resolve):")
for name, ff in tidal_freqs.items():
    # Check if PSD has any power at this frequency
    idx = np.argmin(np.abs(f - ff))
    print(f"    {name:20s} = {ff:.2e} Hz → PSD bin at {f[idx]:.2e} Hz: "
          f"P424={psd_424[idx]:.3e}, P464={psd_464[idx]:.3e}")

# Compare LF spectral slope pre vs post
p424_post_1hz = sp_signal.decimate(p424[post_mask], dec_factor, ftype='fir') if np.sum(post_mask) > dec_factor*60 else None

if p424_post_1hz is not None and len(p424_post_1hz) > 600:
    f_post, psd_424_post = sp_signal.welch(p424_post_1hz, fs=1, nperseg=600, noverlap=300)
    
    # Spectral slope in 0.001-0.1 Hz band
    lf_mask = (f >= 0.001) & (f <= 0.1)
    lf_mask_post = (f_post >= 0.001) & (f_post <= 0.1)
    
    if np.sum(lf_mask) > 5 and np.sum(lf_mask_post) > 5:
        # Fit log-log slope
        log_f = np.log10(f[lf_mask])
        log_psd = np.log10(psd_424[lf_mask])
        slope_pre = np.polyfit(log_f, log_psd, 1)[0]
        
        log_f_post = np.log10(f_post[lf_mask_post])
        log_psd_post = np.log10(psd_424_post[lf_mask_post])
        slope_post = np.polyfit(log_f_post, log_psd_post, 1)[0]
        
        print(f"\n  P424 LF spectral slope (0.001-0.1 Hz):")
        print(f"    Pre-seismic:  {slope_pre:.2f}")
        print(f"    Post-seismic: {slope_post:.2f}")
        print(f"    Δslope:       {slope_post - slope_pre:+.2f}")

# =====================================================================
# A4: SPECTROGRAM — Time-frequency analysis
# =====================================================================
print("\n" + "="*60)
print("A4: SPECTROGRAM — Time-frequency structure")
print("="*60)

# For spectrogram, decimate to reduce computation
# Use bandpass filtered 0.5-50 Hz range on full data
p424_bp = bandpass_filter(p424, 0.5, 50, fs)
p464_bp = bandpass_filter(p464, 0.5, 50, fs)

# Decimate to 100 Hz
p424_d = sp_signal.decimate(p424_bp, 2, ftype='fir')
p464_d = sp_signal.decimate(p464_bp, 2, ftype='fir')
fs_d = 100.0
times_d = times[::2]

# Spectrogram parameters
nperseg = int(30 * fs_d)   # 30 second windows
noverlap = int(28 * fs_d)  # 28 second overlap (93% overlap for smoothness)

print(f"  Spectrogram params: window={nperseg/fs_d:.0f}s, overlap={noverlap/fs_d:.0f}s")

# Compute spectrogram for pre-seismic only
pre_mask_d = times_d < EQ_TIME
f_sg, t_sg, Sxx_424 = sp_signal.spectrogram(
    p424_d[pre_mask_d], fs=fs_d, 
    nperseg=nperseg, noverlap=noverlap,
    window='hann', scaling='density'
)

# Find peak frequencies in each time window
print(f"  Spectrogram shape: {Sxx_424.shape} ({len(t_sg)} time points × {len(f_sg)} freq bins)")

# Look for pre-seismic spectral changes
# Divide into early (T-28 to T-20) and late (T-20 to T-0) pre-seismic
t_sg_rel = t_sg + times[0]  # absolute times

early_mask = t_sg_rel < EQ_TIME - 10*60
late_mask = t_sg_rel >= EQ_TIME - 10*60

if np.sum(early_mask) > 0 and np.sum(late_mask) > 0:
    early_psd = np.mean(Sxx_424[:, early_mask], axis=1)
    late_psd = np.mean(Sxx_424[:, late_mask], axis=1)
    
    # Frequency bands to check
    freq_bands_sg = [
        ('ULF', 0.5, 1.0),
        ('VLF', 1.0, 5.0),
        ('LF', 5.0, 15.0),
        ('MF', 15.0, 30.0),
        ('HF', 30.0, 50.0),
    ]
    
    print("\n  Pre-seismic spectral change (last 10 min vs earlier):")
    for name, lo, hi in freq_bands_sg:
        fb_mask = (f_sg >= lo) & (f_sg <= hi)
        if np.sum(fb_mask) > 0:
            early_power = np.mean(early_psd[fb_mask])
            late_power = np.mean(late_psd[fb_mask])
            ratio = late_power / early_power if early_power > 0 else 0
            status = "✅" if abs(ratio - 1) > 0.1 else "—"
            print(f"    {name:5s} [{lo:4.1f}-{hi:4.1f} Hz]: "
                  f"early={early_power:.2e}, late={late_power:.2e}, "
                  f"ratio={ratio:.3f} {status}")

# Also check whole frequency range for tremor signatures
# Look for narrowband spectral peaks that appear only in late pre-seismic
print("\n  Searching for emerging spectral features...")
n_freq = len(f_sg)
emerging_features = []
for fi in range(n_freq):
    early_p = np.mean(Sxx_424[fi, early_mask]) if np.sum(early_mask) > 0 else 0
    late_p = np.mean(Sxx_424[fi, late_mask]) if np.sum(late_mask) > 0 else 0
    if early_p > 0 and late_p / early_p > 2.0:
        emerging_features.append((f_sg[fi], late_p / early_p))

if emerging_features:
    # Sort by ratio descending
    emerging_features.sort(key=lambda x: x[1], reverse=True)
    print(f"  Found {len(emerging_features)} frequency bins with >2× enhancement:")
    for freq, ratio in emerging_features[:10]:
        print(f"    {freq:.2f} Hz: {ratio:.1f}×")
else:
    print("  No significant emerging narrowband features detected")

# =====================================================================
# A3b: CHECK FOR PRE-SEISMIC SPECTRAL ANOMALY AT SPECIFIC TIMES
# Use fine-grained time windows
# =====================================================================
print("\n" + "="*60)
print("A3b: TIME-RESOLVED SPECTRAL EVOLUTION (4 min windows)")
print("="*60)

# 4 min windows, 30s step
win_4m = int(4 * 60 * fs_d)
step_30s = int(30 * fs_d)

# Compute spectrogram segments
n_windows = (len(p424_d[pre_mask_d]) - win_4m) // step_30s + 1

spectral_evolution = []
for i in range(0, len(p424_d[pre_mask_d]) - win_4m, step_30s):
    seg = p424_d[pre_mask_d][i:i+win_4m]
    # Welch PSD for this segment
    f_4m, psd_4m = sp_signal.welch(seg, fs=fs_d, nperseg=min(120, len(seg)//4))
    t_center = times_d[pre_mask_d][i + win_4m//2]
    
    # Power in key bands
    bands_power = {}
    for name, lo, hi in [('ULF', 0.5, 1.0), ('VLF', 1.0, 5.0), ('LF', 5.0, 15.0),
                          ('MF', 15.0, 30.0), ('HF', 30.0, 50.0)]:
        mask = (f_4m >= lo) & (f_4m <= hi)
        if np.sum(mask) > 0:
            bands_power[name] = float(np.mean(psd_4m[mask]))
    
    spectral_evolution.append({
        't_center': float(t_center),
        't_rel': float(t_center - EQ_TIME),
        'bands': bands_power,
        'total_power': float(np.mean(psd_4m))
    })

# Compare last 8 min vs earlier for power change
last_8min = [s for s in spectral_evolution if s['t_center'] >= EQ_TIME - 8*60]
early_20 = [s for s in spectral_evolution if s['t_center'] < EQ_TIME - 15*60]

if last_8min and early_20:
    print("  Spectral power change: last 8 min vs early baseline")
    for band in ['ULF', 'VLF', 'LF', 'MF', 'HF']:
        early_p = np.mean([s['bands'].get(band, 0) for s in early_20])
        late_p = np.mean([s['bands'].get(band, 0) for s in last_8min])
        if early_p > 0 and late_p > 0:
            print(f"    {band:5s}: early={early_p:.2e}, late={late_p:.4e}, "
                  f"ratio={late_p/early_p:.3f}")

# =====================================================================
# SAVE RESULTS
# =====================================================================
print("\n" + "="*60)
print("SAVING RESULTS")
print("="*60)

output = {
    'script': 'ohat_200hz_dirA.py',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()),
    'data': '200 Hz P_f Milun Fault (dn4338tv9g)',
    'eq_time': EQ_TIME,
    'results': {
        'A1_bandpass': results_a1,
        'A2_xcorr': {
            'pre_mean_r': round(pre_r, 4) if pre_xc else None,
            'post_mean_r': round(post_r, 4) if post_xc else None,
            'delta_r': round(post_r - pre_r, 4) if pre_xc and post_xc else None,
        },
        'A4_spectral_evolution': {
            'total_windows': len(spectral_evolution),
            'last_8min_vs_early': {}
        }
    }
}

# Add spectral details
if last_8min and early_20:
    for band in ['ULF', 'VLF', 'LF', 'MF', 'HF']:
        early_p = np.mean([s['bands'].get(band, 0) for s in early_20])
        late_p = np.mean([s['bands'].get(band, 0) for s in last_8min])
        if early_p > 0:
            output['results']['A4_spectral_evolution']['last_8min_vs_early'][band] = {
                'early': round(early_p, 6),
                'late': round(late_p, 6),
                'ratio': round(float(late_p/early_p), 4)
            }

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        elif isinstance(obj, (np.floating,)): return float(obj)
        elif isinstance(obj, (np.ndarray,)): return obj.tolist()
        return super().default(obj)

with open(os.path.join(OUTPUT_DIR, 'ohat_200hz_dirA_results.json'), 'w') as f:
    json.dump(output, f, indent=2, cls=NumpyEncoder)

print(f"Results saved to {OUTPUT_DIR}/ohat_200hz_dirA_results.json")
print(f"Total time: {time.time()-t0:.1f}s")
