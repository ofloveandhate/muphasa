#!/usr/bin/env python
"""Regenerate the snapshot JSONs for the complex-fixture characterisation tests.

Run this when behaviour deliberately changes (e.g. you intentionally alter a
discretisation that shifts grade values). The snapshots live under
``test/snapshots/`` and are exact-value captures, sorted into a canonical
order so file diffs are stable.

Usage:
    python test/update_complex_snapshots.py [fixture_name ...]
"""
import json
import sys
from pathlib import Path

import mph
from _complex_fixtures import COMPLEX_FIXTURES

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"


def _canonical_pairings(pairings):
    return sorted(
        ([list(birth), sorted(list(s) for s in syzygies)] for birth, syzygies in pairings),
        key=lambda p: (p[0], p[1]),
    )


def regenerate(name):
    builder = COMPLEX_FIXTURES[name]
    trajectories, max_metric_value, hom_dim = builder()
    res = mph.compute_spatiotemporal_landscapes_sparse(trajectories, max_metric_value, hom_dim)
    payload = {
        "fixture": name,
        "max_metric_value": max_metric_value,
        "hom_dim": hom_dim,
        "trajectory_shape": list(trajectories.shape),
        "pairings": _canonical_pairings(res.pairings),
    }
    SNAPSHOT_DIR.mkdir(exist_ok=True)
    out_path = SNAPSHOT_DIR / f"{name}.json"
    out_path.write_text(json.dumps(payload, indent=2))
    print(f"  {name}: {len(payload['pairings'])} pairings → {out_path.relative_to(Path.cwd())}")


def main(argv):
    names = argv[1:] or list(COMPLEX_FIXTURES)
    unknown = [n for n in names if n not in COMPLEX_FIXTURES]
    if unknown:
        raise SystemExit(f"unknown fixture(s): {unknown}. Known: {list(COMPLEX_FIXTURES)}")
    for name in names:
        regenerate(name)


if __name__ == "__main__":
    main(sys.argv)
