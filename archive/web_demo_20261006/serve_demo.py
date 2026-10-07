"""Serve the prebuilt demo using only Python's standard library."""
import argparse
import functools
import http.server
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=8080)
    parser.add_argument('--host', default='127.0.0.1')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent / 'demo'
    if not (root / 'index.html').is_file():
        parser.error('Prebuilt demo missing. Run npm ci, npm run build, and python scripts/package_demo.py first.')
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(root))
    with http.server.ThreadingHTTPServer((args.host, args.port), handler) as server:
        print(f'MHW2Coral demo: http://{args.host}:{args.port}', flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass

if __name__ == '__main__':
    main()
