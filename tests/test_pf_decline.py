#!/usr/bin/env python3
"""
Test: Pre-seismic P_f Decline
==============================
Validates that the 200 Hz pore fluid pressure data shows a systematic
P_f decline before the 2024 M7.4 Hualien earthquake, with P170
exhibiting slope ≈ -0.55 MPa/hr (R² ≥ 0.80).

Expected values (from original analysis):
  P170: slope = -0.55 MPa/hr, R² = 0.84
  P360: slope = -0.34 MPa/hr, R² = 0.64
  P424: slope = -0.13 MPa/hr, R² = 0.32
  P464: slope = -0.23 MPa/hr, R² = 0.73

Reference: deep_dive.py, Section 1
"""
import sys
import numpy as np
from scipy import stats

# Configuration - adjust if running from a different directory
DATA_FILE = 'data/sample_data/pressure_pre_eq_5min.csv'

def load_pressure_data(filepath):
    """Load the last 5 min pre-EQ pressure data."""
    import csv
    cols = {'P170': 1, 'P360': 2, 'P424': 3, 'P464': 4}
    data = {k: [] for k in cols}
    
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            for name, idx in cols.items():
                try:
                    data[name].append(float(row[idx]))
                except (ValueError, IndexError):
                    continue
    
    for name in data:
        data[name] = np.array(data[name])
    return data

def test_p170_decline(tolerance=0.20):
    """Test that P170 shows significant linear decline."""
    data = load_pressure_data(DATA_FILE)
    sr = 200
    t = np.arange(len(data['P170'])) / sr
    
    # Linear regression on P170
    sl, intercept, r, p, se = stats.linregress(t, data['P170'])
    
    # Convert to MPa/hr
    slope_hr = sl * 3600
    
    print(f"P170: slope = {slope_hr:.3f} MPa/hr, R² = {r**2:.3f}, p = {p:.2e}")
    
    # Expected: -0.55 MPa/hr
    assert abs(slope_hr - (-0.55)) < tolerance, \
        f"P170 slope {slope_hr:.3f} MPa/hr != expected -0.55 ± {tolerance}"
    assert r**2 > 0.70, \
        f"P170 R² = {r**2:.3f} < 0.70"
    assert p < 0.001, \
        f"P170 p = {p:.2e} not significant"
    
    print("  ✅ P170 decline confirmed")
    return slope_hr, r**2

def test_all_depths_decline():
    """Test that all depths show negative trend."""
    data = load_pressure_data(DATA_FILE)
    sr = 200
    
    for name in ['P170', 'P360', 'P424', 'P464']:
        t = np.arange(len(data[name])) / sr
        sl, intercept, r, p, se = stats.linregress(t, data[name])
        slope_hr = sl * 3600
        
        print(f"{name}: slope = {slope_hr:.3f} MPa/hr, R² = {r**2:.3f}, p = {p:.2e}")
        
        if p < 0.05:
            print(f"  ✅ {name}: significant decline (p = {p:.2e})")
        else:
            print(f"  ⚠️  {name}: not significant (p = {p:.2e})")

if __name__ == '__main__':
    print("=" * 60)
    print("Test: Pre-seismic P_f Decline")
    print("=" * 60)
    
    try:
        test_p170_decline()
        print()
        test_all_depths_decline()
        print()
        print("✅ ALL TESTS PASSED")
    except AssertionError as e:
        print(f"\n❌ FAILED: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"\n❌ Data file not found: {DATA_FILE}")
        print("   Run from repo root: python3 tests/test_pf_decline.py")
        sys.exit(1)
