#!/usr/bin/env python3
import json, urllib.request

# Figshare API
url = "https://api.figshare.com/v2/articles/32895287"
try:
    req = urllib.request.urlopen(url)
    data = json.loads(req.read())
    print(f"Title: {data.get('title','?')}")
    print(f"Published: {data.get('published_date','?')}")
    print(f"Files count: {len(data.get('files',[]))}")
    for f in data.get('files',[]):
        print(f"  {f.get('name','?')} | size={f.get('size','?')} | download_url={f.get('download_url','?')} | link_only={f.get('is_link_only','?')}")
    print(f"\nStatus: {data.get('status','?')}")
    print(f"Is_public: {data.get('is_public','?')}")
    
    # Also check version 4 specifically
    url2 = "https://api.figshare.com/v2/articles/32895287/versions/4"
    req2 = urllib.request.urlopen(url2)
    data2 = json.loads(req2.read())
    print(f"\nVersion 4 files: {len(data2.get('files',[]))}")
    for f in data2.get('files',[]):
        print(f"  {f.get('name','?')} | size={f.get('size','?')} | download_url={f.get('download_url','?')}")
except Exception as e:
    print(f"Error: {e}")
