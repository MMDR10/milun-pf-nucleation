#!/usr/bin/env python3
import json, urllib.request

# Check Mendeley API for file list
url = "https://data.mendeley.com/api/datasets/dn4338tv9g"
try:
    req = urllib.request.urlopen(url)
    data = json.loads(req.read())
    if 'files' in data:
        for f in data['files']:
            print(f.get('name','?'), '|', f.get('size','?'), 'bytes')
    else:
        print("No 'files' key")
        print(str(data)[:2000])
except Exception as e:
    print(f"Error: {e}")

# Also check the specific file listing endpoint
print("\n--- Checking version 2 ---")
url2 = "https://data.mendeley.com/api/datasets/dn4338tv9g/2"
try:
    req2 = urllib.request.urlopen(url2)
    data2 = json.loads(req2.read())
    if 'files' in data2:
        for f in data2['files']:
            print(f.get('name','?'), '|', f.get('size','?'), 'bytes')
    else:
        print(str(data2)[:2000])
except Exception as e:
    print(f"Error: {e}")

print("\n--- Checking public file listing ---")
url3 = "https://data.mendeley.com/public-api/datasets/dn4338tv9g"
try:
    req3 = urllib.request.urlopen(url3)
    data3 = json.loads(req3.read())
    print(str(data3)[:3000])
except Exception as e:
    print(f"Error: {e}")
