#!/usr/bin/env python3
import csv

with open('/tmp/pf_data.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    rows = list(reader)

# Parse time format '7:30:00 AM'
def parse_time(t):
    parts = t.strip().split()
    if not parts:
        return 0
    hms = parts[0].split(':')
    h = int(hms[0])
    m = int(hms[1])
    s = int(hms[2]) if len(hms) > 2 else 0
    if len(parts) > 1 and parts[1] == 'PM' and h < 12:
        h += 12
    if len(parts) > 1 and parts[1] == 'AM' and h == 12:
        h = 0
    return h * 3600 + m * 60 + s

# M7.4 at 07:58:11 local time
eq_time = parse_time('7:58:00 AM')
pre_rows = [r for r in rows if parse_time(r[0]) < eq_time and parse_time(r[0]) > 0]
post_rows = [r for r in rows if parse_time(r[0]) >= eq_time]

def avg(col):
    vals = [float(v) for v in col if v.strip()]
    return sum(vals)/len(vals) if vals else 0

print(f'Pre-EQ rows: {len(pre_rows)}, Post-EQ rows: {len(post_rows)}')

for name, col_idx in [('170m', 1), ('360m', 2), ('424m', 3), ('464m', 4)]:
    pre_avg = avg([r[col_idx] for r in pre_rows])
    post_avg = avg([r[col_idx] for r in post_rows])
    drop = (pre_avg - post_avg) * 100  # MPa to bar
    pct = ((pre_avg - post_avg) / pre_avg) * 100
    print(f'\n=== P_f {name} ===')
    print(f'  Pre-EQ:  {pre_avg:.4f} MPa ({pre_avg*100:.2f} bar)')
    print(f'  Post-EQ: {post_avg:.4f} MPa ({post_avg*100:.2f} bar)')
    print(f'  Drop:    {drop:.2f} bar ({pct:.2f}%)')

# Also get min/max during coseismic period
print('\n=== Coseismic extremes (07:58-08:02) ===')
coseismic = [r for r in rows if 7*3600+58*60 <= parse_time(r[0]) <= 7*3600+62*60]
for name, col_idx in [('170m', 1), ('360m', 2), ('424m', 3), ('464m', 4)]:
    vals = [float(r[col_idx]) for r in coseismic if r[col_idx].strip()]
    print(f'  {name}: min={min(vals)*100:.2f} bar, max={max(vals)*100:.2f} bar, range={(max(vals)-min(vals))*100:.2f} bar')
