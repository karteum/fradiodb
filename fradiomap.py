#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 12:39:58 2026
@author: Adrien Demarez
"""

from http.server import HTTPServer, BaseHTTPRequestHandler
import re
import fradiodb as fr
DBFILE="anfr2026.gpkg"

class MyHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/index.html" or self.path == "/anfr.fgb":
            try:
                with open(self.path, "rb") as f:
                    content = f.read()
                self.send_response(200)
                self.send_header("Content-Type", "text/html" if self.path == "/index.html"  else "application/x-flatgeobuf")
                self.end_headers()
                self.wfile.write(content)
            except FileNotFoundError:
                self.send_error(404)
            return

        # ---- dynamic route /site/<site_id> ----
        match = re.match(r"^/site/([^/]+)$", self.path)
        if match:
            site_id = match.group(1)        
            payload = fr.query(DBFILE, int(site_id), "site")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(payload.encode("utf-8"))
            return

        # ---- otherwise 404 ----
        self.send_error(404)


if __name__ == "__main__":
    server = HTTPServer(("localhost", 8001), MyHandler)
    print("Server running on http://localhost:8001")
    server.serve_forever()