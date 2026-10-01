#!/usr/bin/env python3
"""Build docs/basin-atlas.html from results/basinscan-*.json."""
import json, pathlib

root = pathlib.Path(__file__).resolve().parent.parent
GLYPHS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
KIND_ORDER = {"cycle": 0, "helix": 1, "spiral": 2, "ray": 3, "fixed": 4, "undecided": 5}


def compact(m):
    """Re-index attractors by basin size, merge inert fixed points, cap at 61 classes + other."""
    cells = m["cells"]
    size = {}
    for c in cells:
        size[c] = size.get(c, 0) + 1
    attrs = m["attractors"]
    inert = [i for i, a in enumerate(attrs) if a["kind"] == "fixed" and size[i] == 1]
    real = [i for i in range(len(attrs)) if i not in set(inert)]
    real.sort(key=lambda i: (-size[i], KIND_ORDER[attrs[i]["kind"]], attrs[i]["id"]))
    classes = []
    remap = {}
    for i in real[:60]:
        remap[i] = len(classes)
        classes.append({"kind": attrs[i]["kind"], "id": attrs[i]["id"], "n": size[i]})
    if len(real) > 60:
        other = len(classes)
        classes.append({"kind": "other", "id": [], "n": sum(size[i] for i in real[60:])})
        for i in real[60:]:
            remap[i] = other
    if inert:
        k = len(classes)
        classes.append({"kind": "inert", "id": [], "n": len(inert)})
        for i in inert:
            remap[i] = k
    return {
        "index": m["index"],
        "f": m["features"],
        "classes": classes,
        "cells": "".join(GLYPHS[remap[c]] for c in cells),
    }


def load(path, keep_all):
    d = json.load(open(root / path))
    maps = [compact(m) for m in d["maps"]]
    summary = {
        "universes": len(maps),
        "radius": d["radius"],
        "conic": sum(1 for m in maps if m["f"]["conic"] >= 0.999),
        "single": sum(1 for m in maps if m["f"]["moving"] + m["f"]["sinks"] <= 1 and m["f"]["inert"] <= 1),
        "checker": sum(1 for m in maps if m["f"]["stripes"] >= 0.25),
        "multi": sum(1 for m in maps if m["f"]["moving"] + m["f"]["sinks"] >= 2),
    }
    if not keep_all:
        chosen = {}
        for key, rev in (("moving", True), ("sinks", True), ("boundary", True), ("stripes", True), ("conic", False)):
            ranked = sorted(maps, key=lambda m: (m["f"][key], m["index"]), reverse=rev)
            for m in ranked[:80]:
                chosen[m["index"]] = m
        maps = sorted(chosen.values(), key=lambda m: m["index"])
    return {"summary": summary, "maps": maps, "sampled": not keep_all}


data = {
    "2": load("results/basinscan-size2.json", True),
    "3": load("results/basinscan-size3-sampled.json", False),
}
tmpl = open(root / "scripts/basin-atlas.template.html").read()
# Data goes in its own script file: a multi-megabyte inline <script> sits in the DOM and makes
# every later DOM change slow in browsers that snapshot the document.
blob = "window.BASIN_DATA=" + json.dumps(data, separators=(",", ":")) + ";"
open(root / "docs/basin-atlas-data.js", "w").write(blob)
open(root / "docs/basin-atlas.html", "w").write(tmpl)
print("wrote docs/basin-atlas.html", len(tmpl), "docs/basin-atlas-data.js", len(blob), {k: len(v["maps"]) for k, v in data.items()})
