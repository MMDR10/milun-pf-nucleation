"""
Ô-HAT θ Analysis: Milun Fault P_f pre-seismic anomaly
=====================================================
Compute θ(TopSV, PhaseLock) for the P_f timeseries
"""
import numpy as np
import pandas as pd
from scipy import signal, stats

print("=" * 70)
print("Ô-HAT θ ANALYSIS - Milun Fault 200Hz P_f (2024/04/03 M7.4)")
print("=" * 70)

base = '/app/working/workspaces/tygtDc/data/Milun_Pf/The 200 Hz Sampling Interval of Porewater Pressure'

# ── Load ──
df_p = pd.read_csv(f'{base}/Groundwater Pressure Data .csv')
df_p.columns = ['Time', 'P170', 'P360', 'P424', 'P464']

N = len(df_p)
sr = 200  # Hz
t_arr = np.arange(N) / sr  # seconds from 7:30
eq_time = 28 * 60  # 1680 sec = 7:58 AM

# ── Helper: Ô-HAT θ(TopSV) ──
def theta_topsv(x, y, win_sec=60, sr=200):
    """Compute sliding window Topological SVD angle θ"""
    win = int(win_sec * sr)
    step = int(10 * sr)  # 10-sec step
    
    n_windows = (len(x) - win) // step
    if n_windows < 1:
        return None, None
    
    thetas = []
    times = []
    
    for i in range(n_windows):
        s = i * step
        e = s + win
        xs = x[s:e]
        ys = y[s:e]
        
        # Build delay embedding (τ=50 samples = 0.25s at 200Hz)
        tau = 50
        dim = 5
        def embed(z):
            n = len(z) - tau * (dim - 1)
            return np.column_stack([z[i:n+i] for i in range(0, dim * tau, tau)])
        
        X = embed(xs)
        Y = embed(ys)
        
        if len(X) < 3 or len(Y) < 3:
            continue
        
        # SVD Topological overlap
        Ux, Sx, Vtx = np.linalg.svd(X, full_matrices=False)
        Uy, Sy, Vty = np.linalg.svd(Y, full_matrices=False)
        
        k = min(3, len(Sx), len(Sy))
        Vx_k = Vtx[:k].T
        Vy_k = Vty[:k].T
        
        # Subspace angle via principal angles
        M = Vx_k.T @ Vy_k
        _, s_principal, _ = np.linalg.svd(M)
        # θ = arccos(geometric mean of principal cosines)
        cos_s = np.clip(s_principal, -1, 1)
        theta = np.arccos(np.mean(cos_s[:k]))
        
        thetas.append(theta)
        times.append(s / sr)
    
    return np.array(times), np.array(thetas)

# ── 1. Pre-seismic P_f trend: detrended anomaly ──
pre = df_p.iloc[:int(eq_time * sr)].copy()
print(f"\n📈 PRE-SEISMIC DETRENDED ANALYSIS ({eq_time//60} min)")
print(f"{'Depth':>8s} {'Slope(MPa/s)':>14s} {'R²':>10s} {'Total_drop(MPa)':>16s} {'Total_drop(%)':>14s}")
print("-" * 62)

for d in ['P170', 'P360', 'P424', 'P464']:
    t_pre = t_arr[:len(pre)]
    y = pre[d].values
    slope, intercept, r, p, std = stats.linregress(t_pre, y)
    total_drop = slope * t_pre[-1]
    pct = total_drop / y.mean() * 100
    print(f"{d:>8s} {slope:>14.2e} {r**2:>10.6f} {total_drop:>16.6f} {pct:>14.6f}")

# ── 2. Ô-HAT θ between adjacent depth pairs ──
print(f"\n🔄 Ô-HAT θ(TopSV) — Pre-seismic (7:30-7:58)")
depths = ['P170', 'P360', 'P424', 'P464']
pairs = [('P170','P360'),('P360','P424'),('P424','P464'),('P170','P464')]
win_sizes = [60, 120, 300]  # 1min, 2min, 5min

for w in win_sizes:
    print(f"\n  Window = {w}s")
    print(f"  {'Pair':>12s} {'θ_mean(rad)':>14s} {'θ_std':>10s} {'θ_min':>10s} {'θ_max':>10s}")
    print(f"  {'-'*48}")
    for p1, p2 in pairs:
        pre1 = df_p[p1][:int(eq_time * sr)].values
        pre2 = df_p[p2][:int(eq_time * sr)].values
        times, thetas = theta_topsv(pre1, pre2, win_sec=w)
        if times is not None:
            print(f"  {p1}×{p2:>6s} {thetas.mean():>14.4f} {thetas.std():>10.4f} {thetas.min():>10.4f} {thetas.max():>10.4f}")

# ── 3. Ô-HAT θ trend: is θ increasing before EQ? ──
print(f"\n📉 Ô-HAT θ TREND before earthquake (5-min windows)")
print(f"  Looking for θ increase → phase space divergence signal")

for p1, p2 in pairs:
    pre1 = df_p[p1][:int(eq_time * sr)].values
    pre2 = df_p[p2][:int(eq_time * sr)].values
    times, thetas = theta_topsv(pre1, pre2, win_sec=60)
    if times is not None:
        # Split early vs late
        mid = len(thetas) // 2
        early = thetas[:mid].mean()
        late = thetas[mid:].mean()
        t_stat, p_val = stats.ttest_ind(thetas[:mid], thetas[mid:])
        delta = (late - early) / early * 100
        print(f"  {p1}×{p2:>6s}: early_θ={early:.4f} late_θ={late:.4f} Δ={delta:+.1f}% t={t_stat:.2f} p={p_val:.4f}", end="")
        if p_val < 0.05 and late > early:
            print(" ← 🔥 θ INCREASE!")
        elif p_val < 0.05 and late < early:
            print(" ← 📉 θ decrease")
        else:
            print()

# ── 4. Co-seismic P_f step analysis ──
print(f"\n⚡ CO-SEISMIC STEP (at t={eq_time}s = 07:58:00)")
print(f"  Using 30s pre-EQ vs 30s post-EQ")
print(f"{'Depth':>8s} {'Pre_30s_mean':>14s} {'Post_30s_mean':>14s} {'Step':>12s} {'Z_step':>10s} {'Dir':>8s}")
print("-" * 66)

pre_30s_s = int(eq_time * sr) - 30 * sr
post_30s_s = int(eq_time * sr)
post_30s_e = int(eq_time * sr) + 30 * sr

for d in depths:
    pre_30 = df_p[d].iloc[pre_30s_s:post_30s_s].values
    post_30 = df_p[d].iloc[post_30s_s:post_30s_e].values
    step = post_30.mean() - pre_30.mean()
    z = step / pre_30.std()
    dir_s = "⬆️" if step > 0 else "⬇️"
    print(f"{d:>8s} {pre_30.mean():>14.6f} {post_30.mean():>14.6f} {step:>12.6f} {z:>10.2f} {dir_s:>8s}")

# ── 5. Full timeseries co-seismic analysis ──
print(f"\n🌊 FULL 1-HOUR P_f TREND (subsampled = 1 Hz)")
t_1hz = np.arange(0, 3600)
df_1hz = df_p.iloc[::200].reset_index(drop=True)
print(f"{'Depth':>8s} {'Mean':>12s} {'Trend(MPa/h)':>14s} {'Pre_drop':>12s} {'Post_jump':>12s}")
print("-" * 58)

for d in depths:
    y_1hz = df_1hz[d].values[:3600]
    # Pre seismic
    pre_y = y_1hz[:1680]
    t_pre = np.arange(1680)
    sl_pre, _, _, _, _ = stats.linregress(t_pre, pre_y)
    
    # Post seismic
    post_y = y_1hz[1680:3600]
    t_post = np.arange(len(post_y))
    sl_post, _, _, _, _ = stats.linregress(t_post, post_y)
    
    # Pre-drop amplitude (last 5 min - first 5 min of pre-seismic)
    pre_drop = pre_y[-300] - pre_y[0]
    post_jump = post_y[-1] - post_y[0]
    
    print(f"{d:>8s} {y_1hz.mean():>12.6f} {sl_pre*3600:>14.2e} {pre_drop:>12.6f} {post_jump:>12.6f}")

print(f"\n✅ Deep analysis complete!")
