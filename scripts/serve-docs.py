#!/usr/bin/env python3
"""Serve docs/ on 127.0.0.1:8766 with caching disabled, so rebuilt pages show on a normal reload."""
import functools, http.server, pathlib, sys

class NoCache(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-store, must-revalidate")
        self.send_header("Expires", "0")
        super().end_headers()

port = int(sys.argv[1]) if len(sys.argv) > 1 else 8766
root = pathlib.Path(__file__).resolve().parent.parent / "docs"
handler = functools.partial(NoCache, directory=str(root))
http.server.ThreadingHTTPServer(("127.0.0.1", port), handler).serve_forever()
