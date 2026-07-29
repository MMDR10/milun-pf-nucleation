#!/usr/bin/env python3
import json, urllib.request

# 1. Query CI Taiwan SensorThings API for wells near Milun Fault
# Try the new STA endpoint
base = "http://sta.colife.org.tw/v1.1"

# Search for groundwater related Things
# First, let's check what's at the root
try:
    req = urllib.request.urlopen(f"{base}/Things?$top=5&$select=name,description")
    data = json.loads(req.read())
    print("=== First 5 Things ===")
    for v in data.get('value', []):
        print(f"  {v.get('name','?')}: {str(v.get('description',''))[:80]}")
    print(f"  Total count: {data.get('@iot.count', 'N/A')}")
except Exception as e:
    print(f"Root Error: {e}")

# 2. Try searching for groundwater-specific wells
# The WRA groundwater stations might have "GW" or specific names
search_terms = ['北埔', '新城', '花崗山', '府前', '地下水']
for term in search_terms:
    try:
        url = f"{base}/Things?$filter=contains(name,'{term}')&$top=5"
        req = urllib.request.urlopen(url)
        data = json.loads(req.read())
        vals = data.get('value', [])
        if vals:
            print(f"\n=== '{term}' results ({len(vals)}) ===")
            for v in vals[:5]:
                print(f"  {v.get('name','?')} | id={v.get('@iot.id','?')}")
                # Get the Datastreams
                ds_url = v.get('Datastreams@iot.navigationLink', '')
                if ds_url:
                    try:
                        ds_req = urllib.request.urlopen(ds_url)
                        ds_data = json.loads(ds_req.read())
                        for ds in ds_data.get('value', [])[:2]:
                            print(f"    → Datastream: {ds.get('name','?')} | unit={ds.get('unitOfMeasurement',{}).get('name','?')}")
                            # Get first observation
                            obs_url = ds.get('Observations@iot.navigationLink', '')
                            if obs_url:
                                obs_req = urllib.request.urlopen(f"{obs_url}?$top=2&$orderby=phenomenonTime desc")
                                obs_data = json.loads(obs_req.read())
                                for ob in obs_data.get('value', [])[:2]:
                                    print(f"      → Obs: {ob.get('result','?')} @ {ob.get('phenomenonTime','?')}")
                    except Exception as e2:
                        print(f"    → DS Error: {e2}")
        else:
            print(f"  '{term}': no results")
    except Exception as e:
        print(f"  '{term}' Error: {e}")

# 3. Try to list Datastream types for groundwater
print("\n=== Checking groundwater Datastreams ===")
try:
    url = f"{base}/Datastreams?$top=20&$select=name,unitOfMeasurement"
    req = urllib.request.urlopen(url)
    data = json.loads(req.read())
    for v in data.get('value', []):
        print(f"  {v.get('name','?')} | unit={v.get('unitOfMeasurement',{}).get('name','?')}")
except Exception as e:
    print(f"Error: {e}")
