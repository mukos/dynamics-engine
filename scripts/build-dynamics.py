#!/usr/bin/env python3
"""Build docs/index.html (Dynamics Engine) from the exact census JSON exports (standard start)."""
import json, pathlib
root = pathlib.Path(__file__).resolve().parent.parent
data = {}
for size in (2, 3, 4):
    path = root / (f"results/exact-size{size}.json" if size < 4 else "results/exact-size4-sampled.json")
    try:
        data[str(size)] = json.load(open(path))
    except (OSError, ValueError):
        print(f"skipping size {size}: {path} missing or incomplete")
# The five-fates legend: one size-3 class per fate whose parameter is near a target that animates well.
TARGET = {"fixed": 12, "cycle": 10, "ray": 6, "helix": 6, "spiral": 0}
legend = []
if "3" in data:
    for kind, target in TARGET.items():
        cls = [c for c in data["3"]["classes"] if c["kind"] == kind]
        if cls:
            legend.append(min(cls, key=lambda c: (abs(c["param"] - target), c["index"])))
# Hero rosettes: size-3 universes whose polar plots look good: mid-length cycles, helices, the spiral, a few rays.
hero = []
if "3" in data:
    for c in data["3"]["classes"]:
        k, p = c["kind"], c["param"]
        if (k == "cycle" and 5 <= p <= 40) or k in ("helix", "spiral") or (k == "ray" and 3 <= p <= 12):
            hero.append({key: c[key] for key in ("kind", "param", "index", "pos", "neg")})
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
out = tmpl.replace("__DATA__", json.dumps(data)).replace("__LEGEND__", json.dumps(legend)).replace("__HERO__", json.dumps(hero, separators=(",", ":")))
open(root / "docs/index.html", "w").write(out)
print("wrote docs/index.html", len(out), "sizes", sorted(data), "hero", len(hero))
