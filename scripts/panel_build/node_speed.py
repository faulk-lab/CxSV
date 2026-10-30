#!/usr/bin/env python3
"""Single-thread speed of a node for the panel's Python-heavy callers.

Sniffles2, cuteSV and SVIM are Python (SVIM single-threaded, the long pole of
every sample), so a pure-Python workload is the relevant per-core speed
proxy. Prints: host, best-of-3 seconds for a fixed workload (lower = faster).
Run inside the pipeline image so every node uses the same interpreter:
    cxsv_container <wd> <tmp> python3 scripts/panel_build/node_speed.py
"""
import socket
import time


def work(n=3_000_000):
    d = {}
    acc = 0
    for i in range(n):
        k = (i * 2654435761) & 0xFFFFF
        d[k] = d.get(k, 0) + 1
        acc += k % 7
    s = sorted(d.items())[:1000]
    return acc + len(s)


best = min((lambda t0: (work(), time.perf_counter() - t0)[1])(time.perf_counter()) for _ in range(3))
print(f"{socket.gethostname()}\t{best:.3f}")
