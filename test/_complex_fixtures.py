"""Larger (~1s) trajectory fixtures used by the characterisation suite.

Each entry returns a tuple ``(trajectories, max_metric_value, hom_dim)``.
Trajectory shape is ``(n_agents, n_timesteps, n_dims)`` and dtype is
``float64``. Random fixtures use ``np.random.seed`` (legacy global RNG) so the
output is reproducible across machines.

These fixtures are picked to exercise different topological dynamics — symmetric
rotation, cluster merging, growing variance, brownian noise, regular grids — so
that snapshot regressions point at distinct code paths.
"""
import numpy as np


def swarm_disperse_26x9():
    """26 agents starting tightly clustered, drifting outward over 9 timesteps.
    Closest analogue to the swarm-data use-case the package targets."""
    np.random.seed(1)
    n_agents = 26
    seed_pts = 0.05 * np.random.randn(n_agents, 3)
    trajectories = np.stack(
        [seed_pts + 0.08 * t * np.random.randn(n_agents, 3) for t in range(9)],
        axis=1,
    )
    return trajectories.astype(np.float64), 1.6, 1


def rotating_ring_36x10():
    """36 agents on a noisy ring, rotating through pi/2 over 10 timesteps."""
    np.random.seed(11)
    n_agents = 36
    n_steps = 10
    theta = np.linspace(0, 2 * np.pi, n_agents, endpoint=False) + 0.1 * np.random.randn(n_agents)
    radii = 1.0 + 0.2 * np.random.randn(n_agents)
    z_offset = 0.2 * np.random.randn(n_agents)
    trajectories = np.zeros((n_agents, n_steps, 3))
    for ti, phase in enumerate(np.linspace(0, np.pi / 2, n_steps)):
        trajectories[:, ti, 0] = radii * np.cos(theta + phase)
        trajectories[:, ti, 1] = radii * np.sin(theta + phase)
        trajectories[:, ti, 2] = z_offset * np.cos(phase)
    return trajectories.astype(np.float64), 1.4, 1


def two_clusters_merging_32x9():
    """Two distant clusters of 16 agents each, drifting toward each other and
    eventually overlapping at the final timestep."""
    np.random.seed(2)
    n_per = 16
    n_steps = 9
    cluster_a = 0.20 * np.random.randn(n_per, 3) + np.array([0.0, 0.0, 0.0])
    cluster_b = 0.20 * np.random.randn(n_per, 3) + np.array([1.5, 0.0, 0.0])
    trajectories = np.zeros((2 * n_per, n_steps, 3))
    for ti, alpha in enumerate(np.linspace(0, 1, n_steps)):
        drift = alpha * np.array([0.75, 0.0, 0.0])
        trajectories[:n_per, ti, :] = cluster_a + drift
        trajectories[n_per:, ti, :] = cluster_b - drift
    return trajectories.astype(np.float64), 1.5, 1


def brownian_30x10():
    """30 agents performing independent random walks from random starts."""
    np.random.seed(3)
    n_agents = 30
    n_steps = 10
    increments = 0.15 * np.random.randn(n_agents, n_steps - 1, 3)
    initial = 0.5 * np.random.randn(n_agents, 1, 3)
    trajectories = np.concatenate([initial, initial + np.cumsum(increments, axis=1)], axis=1)
    return trajectories.astype(np.float64), 2.0, 1


def drifting_grid_49x10():
    """7×7 grid in the xy-plane drifting in z by an agent-dependent amount."""
    n_side = 7
    n_steps = 10
    xs, ys = np.meshgrid(np.linspace(0, 1, n_side), np.linspace(0, 1, n_side), indexing="ij")
    positions_xy = np.stack([xs.ravel(), ys.ravel()], axis=-1).astype(np.float64)
    n_agents = n_side * n_side
    trajectories = np.zeros((n_agents, n_steps, 3))
    for ti, dz in enumerate(np.linspace(0, 0.6, n_steps)):
        trajectories[:, ti, 0:2] = positions_xy
        trajectories[:, ti, 2] = dz * np.cos(np.arange(n_agents))
    return trajectories.astype(np.float64), 0.55, 1


COMPLEX_FIXTURES = {
    "swarm_disperse_26x9":       swarm_disperse_26x9,
    "rotating_ring_36x10":       rotating_ring_36x10,
    "two_clusters_merging_32x9": two_clusters_merging_32x9,
    "brownian_30x10":            brownian_30x10,
    "drifting_grid_49x10":       drifting_grid_49x10,
}
