import urllib.request, json, csv, os

base = '/app/working/workspaces/tygtDc/data/Milun_Pf'

# Download groundwater monitoring well data API
url = 'https://data.moenv.gov.tw/api/v2/gisepa_p_33?limit=1000&format=JSON'
try:
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.loads(resp.read().decode())
        print(f"API returned type: {type(data).__name__}")
        print(f"Keys: {list(data.keys()) if isinstance(data, dict) else 'list'}")
        if isinstance(data, dict):
            for k, v in data.items():
                if isinstance(v, list):
                    print(f"  {k}: {len(v)} items")
                    if len(v) > 0:
                        print(f"  First item keys: {list(v[0].keys()) if isinstance(v[0], dict) else type(v[0])}")
                else:
                    print(f"  {k}: {v}")
except Exception as e:
    print(f"Error: {e}")

# Try alternative: the well condition dataset
url2 = 'https://data.gov.tw/dataset/32718'
print(f"\nDataset page: {url2}")

# Try the Water Resources Agency API
print("\n--- Trying WRA groundwater data ---")
wra_urls = [
    'https://gweb.wra.gov.tw/HydroInfoCore/DataIntroduction?data=%E5%9C%B0%E4%B8%8B%E6%B0%B4%E4%BD%8D',
]
for u in wra_urls:
    try:
        req = urllib.request.Request(u, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            content = resp.read().decode('utf-8', errors='replace')[:500]
            print(f"\nWRA data page ({len(content)} chars): {content[:200]}")
    except Exception as e:
        print(f"Error accessing WRA: {e}")
