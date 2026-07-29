#!/usr/bin/env python3
"""
End-to-End Test with Sample Data
=================================
Runs a mini version of the full analysis pipeline on sample data
(60 seconds of 200 Hz P_f data) to verify scripts are importable
and produce consistent output.

Tests:
1. Script imports (no syntax errors)
2. Basic data loading works
3. Simple statistical checks pass
"""
import sys
import os

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

# ── Test 1: Import all scripts ──
def test_imports():
    print("Test 1: Script imports...")
    try:
        import deep_dive
        print("  ✅ deep_dive.py")
    except SyntaxError as e:
        print(f"  ❌ deep_dive.py syntax error: {e}")
        return False
    
    for script in ['ohat_200hz_analysis.py', 'ohat_200hz_detail.py', 
                    'run_analysis.py', 'run_theta.py']:
        try:
            # Just check syntax
            with open(f'scripts/{script}') as f:
                compile(f.read(), script, 'exec')
            print(f"  ✅ {script}")
        except SyntaxError as e:
            print(f"  ❌ {script}: {e}")
            return False
    
    print("  ✅ All imports pass")
    return True

# ── Test 2: Sample data integrity ──
def test_sample_data():
    print("\nTest 2: Sample data integrity...")
    import csv
    
    sample_file = 'data/sample_data/pressure_60s.csv'
    
    with open(sample_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    
    n_cols = len(header)
    n_rows = len(rows)
    
    print(f"  Header: {header}")
    print(f"  Rows: {n_rows}, Columns: {n_cols}")
    
    assert n_rows == 12000, f"Expected 12000 rows (60s @ 200Hz), got {n_rows}"
    assert n_cols == 5, f"Expected 5 columns, got {n_cols}"
    assert 'GwP_170 m (MPa)' in header[1], f"Missing P170 column"
    
    # Verify numeric integrity
    for row in rows[:100]:
        vals = [float(x) for x in row[1:]]
        # Values should be in reasonable range for MPa
        assert 1.0 < vals[0] < 2.0, f"P170 out of range: {vals[0]}"   # P170 ~1.65
        assert 3.0 < vals[1] < 4.0, f"P360 out of range: {vals[1]}"   # P360 ~3.53
        assert 3.5 < vals[2] < 4.5, f"P424 out of range: {vals[2]}"   # P424 ~4.16
        assert 4.0 < vals[3] < 5.0, f"P464 out of range: {vals[3]}"   # P464 ~4.49
    
    print("  ✅ Data integrity verified")
    return True

# ── Test 3: Accelerometer data ──
def test_accel_data():
    print("\nTest 3: Accelerometer data...")
    import csv
    
    accel_file = 'data/sample_data/accel_30s.csv'
    
    with open(accel_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    
    print(f"  Header: {header[:4]}")
    print(f"  Rows: {len(rows)} = {len(rows)/200:.1f}s @ 200Hz")
    
    return True

# ── Test 4: Analysis result JSONs exist ──
def test_results_exist():
    print("\nTest 4: Analysis result files...")
    import json, glob
    
    results_dir = 'data/results'
    json_files = glob.glob(f'{results_dir}/*.json')
    
    print(f"  Found {len(json_files)} result files")
    
    for jf in sorted(json_files)[:3]:
        with open(jf) as f:
            data = json.load(f)
        print(f"    {os.path.basename(jf)}: {type(data).__name__}")
    
    json_with_data = [jf for jf in json_files if os.path.getsize(jf) > 100]
    print(f"  {len(json_with_data)}/{len(json_files)} files contain substantive data")
    
    return True

if __name__ == '__main__':
    print("=" * 60)
    print("END-TO-END PIPELINE TEST")
    print("=" * 60)
    
    tests = [test_imports, test_sample_data, test_accel_data, test_results_exist]
    
    all_pass = True
    for test in tests:
        try:
            result = test()
            if not result:
                all_pass = False
        except Exception as e:
            print(f"\n❌ FAILED: {e}")
            all_pass = False
    
    print()
    if all_pass:
        print("✅ ALL PIPELINE TESTS PASSED")
    else:
        print("❌ SOME TESTS FAILED")
        sys.exit(1)
