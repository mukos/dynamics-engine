#!/usr/bin/env python3
"""Merge JSON outputs of sharded runs (`--json --offset i --stride k`) into one histogram."""
import json, sys

merged, total, meta = {}, 0, None
for path in sys.argv[1:]:
    with open(path) as f:
        data = json.load(f)
    key = (data["size"], data["steps"], data["skip_previous"])
    if meta is None:
        meta = key
    elif meta != key:
        sys.exit(f"{path}: size/steps/rule {key} differs from {meta}")
    for sig, n in data["results"].items():
        merged[sig] = merged.get(sig, 0) + n
    total += data["total"]

rows = sorted(merged.items(), key=lambda kv: (-kv[1], kv[0]))
print(json.dumps({"size": meta[0], "steps": meta[1], "skip_previous": meta[2], "total": total,
                  "shards": len(sys.argv) - 1, "results": dict(rows)}, indent=2))
