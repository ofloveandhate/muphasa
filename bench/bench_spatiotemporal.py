#!/usr/bin/env python
"""Benchmark harness for compute_spatiotemporal_landscapes_sparse.

Each fixture runs in its own subprocess so peak-RSS is isolated (the value
returned by getrusage is a high-water mark for the whole process and never
shrinks; running in-process would only let us see the largest fixture's peak).

Usage:
    python bench/bench_spatiotemporal.py
    python bench/bench_spatiotemporal.py --save bench/baseline.json
    python bench/bench_spatiotemporal.py --compare bench/baseline.json
    python bench/bench_spatiotemporal.py --fixtures rand_20x8 rand_25x10
"""
import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


# ---------- fixtures ----------

def _stationary_unit_square():
    pts = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]])
    return np.repeat(pts[:, None, :], 2, axis=1).astype(np.float64), 2.0, 1


def _random_trajectory(n_agents, n_steps, seed, max_metric_value=1.5, hom_dim=1):
    np.random.seed(seed)
    return np.random.rand(n_agents, n_steps, 3).astype(np.float64), max_metric_value, hom_dim


FIXTURES = {
    "square_4x2":   _stationary_unit_square,
    "rand_10x5":   lambda: _random_trajectory(10, 5,  seed=123),
    "rand_15x8":   lambda: _random_trajectory(15, 8,  seed=123),
    "rand_20x8":   lambda: _random_trajectory(20, 8,  seed=123),
    "rand_25x10":  lambda: _random_trajectory(25, 10, seed=123),
    "rand_30x10":  lambda: _random_trajectory(30, 10, seed=123),
}


# ---------- worker: runs one fixture, prints JSON ----------

def _worker(name):
    import mph
    trajectories, max_metric_value, hom_dim = FIXTURES[name]()
    t0 = time.perf_counter()
    res = mph.compute_spatiotemporal_landscapes_sparse(trajectories, max_metric_value, hom_dim)
    elapsed = time.perf_counter() - t0
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        peak_kb //= 1024  # macOS reports bytes; Linux reports KB

    print(json.dumps({
        "fixture": name,
        "time_s": elapsed,
        "peak_rss_kb": peak_kb,
        "n_pairings": len(res.pairings),
    }))


# ---------- driver ----------

def _run_fixture(name, repeat):
    runs = []
    for _ in range(repeat):
        out = subprocess.check_output([sys.executable, __file__, "_worker", name])
        runs.append(json.loads(out))
    times = sorted(r["time_s"] for r in runs)
    return {
        "time_s_min":    times[0],
        "time_s_median": times[len(times) // 2],
        "time_s_max":    times[-1],
        "peak_rss_kb":   max(r["peak_rss_kb"] for r in runs),
        "n_pairings":    runs[0]["n_pairings"],
    }


def _print_table(results, baseline=None):
    cols = ["fixture", "time(ms)", "min/max(ms)", "peak(MB)", "pairings"]
    if baseline:
        cols += ["Δtime", "Δrss"]
    rows = []
    for name, r in results.items():
        time_ms = r["time_s_median"] * 1000
        time_range = f"{r['time_s_min']*1000:.1f}/{r['time_s_max']*1000:.1f}"
        rss_mb = r["peak_rss_kb"] / 1024
        row = [name, f"{time_ms:.1f}", time_range, f"{rss_mb:.1f}", str(r["n_pairings"])]
        if baseline and name in baseline:
            b = baseline[name]
            dt = (r["time_s_median"] - b["time_s_median"]) / b["time_s_median"] * 100
            dr = (r["peak_rss_kb"] - b["peak_rss_kb"]) / b["peak_rss_kb"] * 100
            row += [f"{dt:+.1f}%", f"{dr:+.1f}%"]
        elif baseline:
            row += ["-", "-"]
        rows.append(row)

    widths = [max(len(c), *(len(r[i]) for r in rows)) for i, c in enumerate(cols)]
    fmt = "  ".join(f"{{:>{w}}}" for w in widths)
    print(fmt.format(*cols))
    print(fmt.format(*("-" * w for w in widths)))
    for row in rows:
        print(fmt.format(*row))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--save", help="write results to this JSON path")
    p.add_argument("--compare", help="diff against a previously-saved JSON")
    p.add_argument("--fixtures", nargs="*", help="fixtures to run (default: all)")
    p.add_argument("--repeat", type=int, default=3, help="iterations per fixture (default: 3)")
    args = p.parse_args()

    names = args.fixtures or list(FIXTURES)
    unknown = [n for n in names if n not in FIXTURES]
    if unknown:
        raise SystemExit(f"unknown fixture(s): {unknown}. Known: {list(FIXTURES)}")

    baseline = json.loads(Path(args.compare).read_text()) if args.compare else None

    results = {}
    for name in names:
        print(f"running {name} ...", file=sys.stderr, flush=True)
        results[name] = _run_fixture(name, args.repeat)

    print()
    _print_table(results, baseline)

    if args.save:
        Path(args.save).write_text(json.dumps(results, indent=2))
        print(f"\nsaved → {args.save}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "_worker":
        _worker(sys.argv[2])
    else:
        main()
