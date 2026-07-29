import pandas as pd

path = '/app/working/workspaces/tygtDc/data/Milun_Pf/The 200 Hz Sampling Interval of Porewater Pressure/Groundwater Pressure Data .csv'

df = pd.read_csv(path, nrows=200000)

unique_times = df['Datetime'].unique()
print(f"Unique datetime values: {len(unique_times)}")
for t in unique_times[:20]:
    print(f"  '{t}'")
if len(unique_times) > 20:
    print(f"  ... and {len(unique_times)-20} more")

print(f"\nRow 0:    {df.iloc[0]['Datetime']} | P_f[170]={df.iloc[0]['GwP_170 m (MPa)']:.6f}")
print(f"Row 1000:  time={df.iloc[1000]['Datetime']}")
print(f"Row 36000: time={df.iloc[36000]['Datetime']}")

eq_row = 28 * 60 * 200
print(f"\nExpected EQ row: ~{eq_row}")

r2 = df.iloc[300000:310000] if 310000 < len(df) else df.iloc[-10000:]
print(f"\nRows ~300k: unique times = {r2['Datetime'].nunique()}")
print(f"First in batch: {r2.iloc[0]['Datetime']}")
print(f"Last in batch: {r2.iloc[-1]['Datetime']}")
