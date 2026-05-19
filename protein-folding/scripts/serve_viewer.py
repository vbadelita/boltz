#!/usr/bin/env python3
from __future__ import annotations

import argparse
import http.server
import socketserver
from pathlib import Path

VIEWER_ROOT = Path(__file__).resolve().parents[1] / "viewer"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Serve the protein-folding viewer.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8766)
    return parser.parse_args()


class ViewerRootHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(VIEWER_ROOT), **kwargs)


class ReusableTCPServer(socketserver.TCPServer):
    allow_reuse_address = True


def main() -> int:
    args = parse_args()
    with ReusableTCPServer((args.host, args.port), ViewerRootHandler) as httpd:
        print(f"Serving {VIEWER_ROOT} at http://{args.host}:{args.port}/")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            return 0


if __name__ == "__main__":
    raise SystemExit(main())
