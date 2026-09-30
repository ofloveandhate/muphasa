"""Ctrl-C must be able to stop a running landscape computation.

The C++ code holds the GIL for the whole call, so Python only gets to raise KeyboardInterrupt
if the C++ loops poll for pending signals. A timer signal stands in for the keypress here.
"""

import signal
import time

import numpy as np
import pytest

from mph import compute_spatiotemporal_landscapes_tree_sparse

pytestmark = pytest.mark.skipif(
    not hasattr(signal, "setitimer"), reason="needs signal.setitimer"
)


def _tree_input(n_agents, n_steps, seed=123):
    """Tree-structured input with no branching: every agent has one lineage."""
    rng = np.random.default_rng(seed)
    positions = [rng.random((n_agents, 3)) for _ in range(n_steps)]
    parents = [np.full(n_agents, -1, dtype=np.int32)] + [
        np.arange(n_agents, dtype=np.int32) for _ in range(n_steps - 1)
    ]
    return positions, parents


def _raise_keyboard_interrupt(signum, frame):
    raise KeyboardInterrupt


def test_tree_sparse_is_interruptible():
    small = _tree_input(8, 4)
    expected = compute_spatiotemporal_landscapes_tree_sparse(*small, 1.5, 1).pairings

    # Several seconds of work if left to run to completion.
    positions, parents = _tree_input(30, 10)
    delay = 0.5

    old_handler = signal.signal(signal.SIGALRM, _raise_keyboard_interrupt)
    try:
        signal.setitimer(signal.ITIMER_REAL, delay)
        start = time.perf_counter()
        with pytest.raises(KeyboardInterrupt):
            compute_spatiotemporal_landscapes_tree_sparse(positions, parents, 1.5, 1)
        elapsed = time.perf_counter() - start
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)

    # Without polling, the exception only surfaces once the whole computation has finished.
    assert elapsed < delay + 1.5

    # An interrupted call must not poison later ones.
    again = compute_spatiotemporal_landscapes_tree_sparse(*small, 1.5, 1).pairings
    assert again == expected
