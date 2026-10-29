"""Parallel ranged download (sandbox) of full summary-statistic files for overlap correlation (round 6)."""
import sys, os, urllib.request, concurrent.futures as cf
url, out, nthr = sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 32
req = urllib.request.Request(url, method="HEAD"); size = int(urllib.request.urlopen(req, timeout=60).headers["Content-Length"])
chunk = -(-size // nthr); parts = [(i, i * chunk, min(size, (i + 1) * chunk) - 1) for i in range(nthr)]
def get(p):
    i, a, b = p; fn = f"{out}.part{i:03d}"
    if os.path.exists(fn) and os.path.getsize(fn) == b - a + 1: return i
    for attempt in range(6):
        try:
            r = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
            with urllib.request.urlopen(r, timeout=300) as resp, open(fn, "wb") as f:
                while True:
                    d = resp.read(1 << 20)
                    if not d: break
                    f.write(d)
            if os.path.getsize(fn) == b - a + 1: return i
        except Exception as e:
            print("retry", i, e, flush=True)
    raise RuntimeError(f"part {i} failed")
with cf.ThreadPoolExecutor(nthr) as ex: list(ex.map(get, parts))
with open(out, "wb") as f:
    for i, a, b in parts:
        fn = f"{out}.part{i:03d}"; f.write(open(fn, "rb").read()); os.remove(fn)
print("done", out, os.path.getsize(out), "expected", size)
