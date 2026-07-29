import pandas as pd
import os

base = '/app/working/workspaces/tygtDc/data/Milun_Pf'

for root, dirs, files in os.walk(base):
    for f in files:
        if f.endswith('.csv'):
            path = os.path.join(root, f)
            df = pd.read_csv(path, nrows=10)
            print(f'=== {f} ===')
            print(f'Path: {path}')
            print(f'Columns: {list(df.columns)}')
            print(df.head(10).to_string())
            print()
            
            # Count total rows
            with open(path) as fh:
                for i,l in enumerate(fh): pass
            print(f'Total rows (excl header): {i}')
            print()
