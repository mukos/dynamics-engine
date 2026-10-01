#!/usr/bin/env python3
"""Build docs/index.html (Universal Dynamics) from the exact census JSON exports (standard start)."""
import json, pathlib
root = pathlib.Path(__file__).resolve().parent.parent
data = {}
for size in (2, 3, 4):
    path = root / (f"results/exact-size{size}.json" if size < 4 else "results/exact-size4-sampled.json")
    try:
        data[str(size)] = json.load(open(path))
    except (OSError, ValueError):
        print(f"skipping size {size}: {path} missing or incomplete")
# Keep only the classes the page shows: the most common and the 8 rarest per fate.
for d in data.values():
    keep = []
    for kind in ("fixed", "cycle", "ray", "helix", "spiral", "undecided"):
        cls = [c for c in d["classes"] if c["kind"] == kind]
        if not cls:
            continue
        cls.sort(key=lambda c: (-c["count"], c["param"]))
        typical, rest = cls[0], cls[1:]
        rest.sort(key=lambda c: (c["count"], -c["param"]))
        keep.append(typical)
        keep.extend(rest[:8])
    d["classes"] = keep
tmpl = open(root / "scripts/dynamics.template.html").read()
out = tmpl.replace("__DATA__", json.dumps(data))
open(root / "docs/index.html", "w").write(out)
print("wrote docs/index.html", len(out), "sizes", sorted(data))
