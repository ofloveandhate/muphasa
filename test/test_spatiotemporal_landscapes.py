"""Characterisation tests for the spatiotemporal landscape pipeline.

These freeze the observable behaviour of `compute_spatiotemporal_landscapes_sparse`
(plus the `SparseLandscape` and `Grade` helpers) so we can refactor the C++
internals without silent regressions. Snapshots are exact-value captures from
the current implementation; the tolerance in `_assert_pairings_equal` is for
float robustness across platforms, not because the expected values are
approximate.
"""

import json
import math
import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

from mph import SparseLandscape, compute_spatiotemporal_landscapes_sparse
from mph._mph import Grade

from _complex_fixtures import COMPLEX_FIXTURES

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"


# ---------- helpers ----------

SQRT2 = math.sqrt(2.0)


def _canonical(pairings):
    """Sort pairings (and syzygies within each) into a canonical order so that
    the comparison doesn't depend on the C++ output ordering."""
    return sorted(
        ((tuple(birth), sorted(tuple(s) for s in syzygies)) for birth, syzygies in pairings),
        key=lambda p: (p[0], p[1]),
    )


def _assert_pairings_equal(actual, expected, tol=1e-9):
    """Compare two lists of (birth, syzygies) up to canonical ordering and
    a small float tolerance."""
    a = _canonical(actual)
    e = _canonical(expected)
    assert len(a) == len(e), f"pairing count mismatch: got {len(a)}, expected {len(e)}\nactual={a}\nexpected={e}"
    for (ab, asy), (eb, esy) in zip(a, e):
        assert ab == pytest.approx(eb, abs=tol), f"birth mismatch: {ab} vs {eb}"
        assert len(asy) == len(esy), f"syzygy count for birth {eb} mismatch"
        for asy_i, esy_i in zip(asy, esy):
            assert asy_i == pytest.approx(esy_i, abs=tol), f"syzygy mismatch: {asy_i} vs {esy_i}"


def _stationary_trajectory(points, n_timesteps):
    """Repeat a (n_agents, n_dims) point cloud across n_timesteps."""
    points = np.asarray(points, dtype=np.float64)
    return np.repeat(points[:, None, :], n_timesteps, axis=1)


# ---------- fixtures (small, fast inputs) ----------

@pytest.fixture
def stationary_unit_square():
    """Four corners of a unit square in 3D, stationary over 2 timesteps.

    Has H_1 = Z because the four points form a 1-cycle that is not filled
    until distance reaches the diagonal sqrt(2).
    """
    return _stationary_trajectory(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        n_timesteps=2,
    )


@pytest.fixture
def stationary_triangle():
    """Three corners of an equilateral triangle, stationary."""
    return _stationary_trajectory(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, math.sqrt(3) / 2, 0.0]],
        n_timesteps=2,
    )


@pytest.fixture
def stationary_tetrahedron():
    """Four corners of a tetrahedron — fills immediately so H_1 is trivial."""
    return _stationary_trajectory(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        n_timesteps=3,
    )


@pytest.fixture
def drifting_square():
    """Square that grows from side=1 at t=0 to side=2 at t=2."""
    def square(side):
        return [[0.0, 0.0, 0.0], [side, 0.0, 0.0], [side, side, 0.0], [0.0, side, 0.0]]
    sides = [1.0, 1.5, 2.0]
    # shape (n_agents=4, n_timesteps=3, n_dims=3)
    return np.stack([np.array(square(s)) for s in sides], axis=1).astype(np.float64)


@pytest.fixture
def random_small_seeded():
    """5 agents, 4 timesteps, 3D — fixed seed via legacy global RNG to keep the
    snapshot reproducible across machines (numpy preserves the legacy
    Mersenne-Twister stream)."""
    np.random.seed(42)
    return np.random.rand(5, 4, 3).astype(np.float64)


# ---------- tests: shape / sanity ----------

def test_returns_sparse_landscape(stationary_unit_square):
    res = compute_spatiotemporal_landscapes_sparse(stationary_unit_square, max_metric_value=2.0, hom_dim=1)
    assert isinstance(res, SparseLandscape)
    # pairings is a list; each element is (Grade-like, [Grade-like, ...])
    for birth, syzygies in res.pairings:
        assert len(birth) == 3
        for s in syzygies:
            assert len(s) == 3


def test_triangle_has_no_H1(stationary_triangle):
    res = compute_spatiotemporal_landscapes_sparse(stationary_triangle, max_metric_value=2.0, hom_dim=1)
    assert res.pairings == []


def test_tetrahedron_has_no_H1(stationary_tetrahedron):
    res = compute_spatiotemporal_landscapes_sparse(stationary_tetrahedron, max_metric_value=2.0, hom_dim=1)
    assert res.pairings == []


# ---------- tests: snapshot pairings ----------

def test_stationary_unit_square_snapshot(stationary_unit_square):
    res = compute_spatiotemporal_landscapes_sparse(stationary_unit_square, max_metric_value=2.0, hom_dim=1)
    expected = [
        ((1.0, 0.0, 1.0), [(1.0, 0.0, SQRT2)]),
        ((0.0, 1.0, 1.0), [(1.0, 1.0, 1.0), (0.0, 1.0, SQRT2)]),
    ]
    _assert_pairings_equal(res.pairings, expected)


def test_drifting_square_snapshot(drifting_square):
    res = compute_spatiotemporal_landscapes_sparse(drifting_square, max_metric_value=3.0, hom_dim=1)
    expected = [
        ((2.0, 0.0, 1.0), [(2.0, 0.0, SQRT2)]),
        ((1.0, 1.0, 1.5), [(2.0, 1.0, 1.5), (1.0, 1.0, 1.5 * SQRT2)]),
        ((0.0, 2.0, 2.0), [(1.0, 2.0, 2.0), (0.0, 2.0, 2.0 * SQRT2)]),
    ]
    _assert_pairings_equal(res.pairings, expected)


def test_random_seeded_5x4x3_snapshot(random_small_seeded):
    res = compute_spatiotemporal_landscapes_sparse(random_small_seeded, max_metric_value=1.5, hom_dim=1)
    expected = [
        (
            (1.0, 2.0, 0.7517775532690216),
            [
                (2.0, 2.0, 0.7517775532690216),
                (1.0, 3.0, 0.7517775532690216),
                (1.0, 2.0, 0.8361581100771492),
            ],
        ),
        (
            (3.0, 0.0, 0.8201421885341216),
            [
                (3.0, 1.0, 0.8201421885341216),
                (3.0, 0.0, 0.8572349627900021),
            ],
        ),
    ]
    _assert_pairings_equal(res.pairings, expected)


def test_pipeline_is_deterministic(stationary_unit_square):
    """Running twice on the same input must give identical output."""
    a = compute_spatiotemporal_landscapes_sparse(stationary_unit_square, max_metric_value=2.0, hom_dim=1)
    b = compute_spatiotemporal_landscapes_sparse(stationary_unit_square, max_metric_value=2.0, hom_dim=1)
    _assert_pairings_equal(a.pairings, b.pairings, tol=0.0)


# ---------- tests: SparseLandscape.eval ----------

@pytest.fixture
def square_landscape(stationary_unit_square):
    return compute_spatiotemporal_landscapes_sparse(stationary_unit_square, max_metric_value=2.0, hom_dim=1)


@pytest.mark.parametrize(
    "grade,k,expected",
    [
        # k <= 0 → 0
        ((1.0, 0.0, 1.0), 0, 0),
        ((1.0, 0.0, 1.0), -1, 0),
        # at the birth grade itself: persistence = 0
        ((1.0, 0.0, 1.0), 1, 0.0),
        # part-way through the bar
        ((1.0, 0.0, 1.2), 1, pytest.approx(0.2, abs=1e-9)),
        # only one bar live at this grade -> layer 2 is 0
        ((1.0, 0.0, 1.2), 2, 0),
        # grade exceeds the syzygy in scale -> bar already dead
        ((1.0, 1.0, 1.5), 1, 0),
        # second bar live with small persistence
        ((0.0, 1.0, 1.4), 1, pytest.approx(SQRT2 - 1.4, abs=1e-9)),
        # grade not <= any birth → 0
        ((-1.0, -1.0, -1.0), 1, 0),
    ],
)
def test_eval_at_grades(square_landscape, grade, k, expected):
    assert square_landscape.eval(grade, k) == expected


def test_eval_above_layer_count_returns_zero(square_landscape):
    # Only at most 2 bars can ever be live at the same grade for this fixture
    assert square_landscape.eval((1.0, 1.0, 5.0), 5) == 0


# ---------- tests: SparseLandscape.eval2 (batched C++ eval) ----------
#
# eval2 is the batched C++ replacement for eval. The per-cell Python eval is
# the spec, so the cross-check tests below assert eval2 equals the per-cell
# loop on every fixture we own. Edge cases and structural invariants follow.

@pytest.fixture
def drifting_square_landscape(drifting_square):
    return compute_spatiotemporal_landscapes_sparse(drifting_square, max_metric_value=3.0, hom_dim=1)


@pytest.fixture
def triangle_landscape(stationary_triangle):
    """Empty pairings — exercises the no-pairings path."""
    return compute_spatiotemporal_landscapes_sparse(stationary_triangle, max_metric_value=2.0, hom_dim=1)


@pytest.fixture
def random_landscape(random_small_seeded):
    return compute_spatiotemporal_landscapes_sparse(random_small_seeded, max_metric_value=1.5, hom_dim=1)


def _eval_cellwise(landscape, grades, k_max):
    """Reference: per-cell Python eval over a (N, d) grid for layers 1..k_max."""
    return np.array(
        [[landscape.eval(tuple(g), k) for k in range(1, k_max + 1)] for g in grades]
    )


def _grid_around_landscape(landscape, n_per_axis=4, pad=0.2):
    """Build a Cartesian grid that brackets the landscape's birth/syzygy support
    on every axis, with a small pad so we test the regions just outside support too."""
    if not landscape.pairings:
        # Arbitrary single grade — the empty-landscape tests pick their own grid anyway.
        return np.zeros((1, 3))
    all_grades = []
    for birth, syzygies in landscape.pairings:
        all_grades.append(list(birth))
        all_grades.extend(list(s) for s in syzygies)
    arr = np.array(all_grades)
    mins, maxs = arr.min(axis=0), arr.max(axis=0)
    axes = [np.linspace(mins[i] - pad, maxs[i] + pad, n_per_axis) for i in range(arr.shape[1])]
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, arr.shape[1])


def test_eval2_matches_eval_on_grid(square_landscape):
    """Cross-check the batched C++ path against the per-cell Python path on a small grid."""
    t_grid = np.linspace(-0.5, 2.5, 7)
    s_grid = np.linspace(0.0, 2.0, 9)
    k_max = 3

    grades = np.stack(np.meshgrid(t_grid, t_grid, s_grid, indexing="ij"), axis=-1).reshape(-1, 3)
    out = square_landscape.eval2(grades, k_max)

    assert out.shape == (grades.shape[0], k_max)
    np.testing.assert_allclose(out, _eval_cellwise(square_landscape, grades, k_max), atol=1e-12)


@pytest.mark.parametrize(
    "landscape_fixture",
    ["square_landscape", "drifting_square_landscape", "random_landscape"],
)
def test_eval2_matches_eval_per_fixture(landscape_fixture, request):
    """Cross-check on each small landscape using a grid auto-sized to its support."""
    landscape = request.getfixturevalue(landscape_fixture)
    grades = _grid_around_landscape(landscape, n_per_axis=5)
    k_max = 3
    out = landscape.eval2(grades, k_max)
    np.testing.assert_allclose(out, _eval_cellwise(landscape, grades, k_max), atol=1e-12)


@pytest.mark.parametrize(
    "grade,k,expected",
    [
        ((1.0, 0.0, 1.0), 1, 0.0),
        ((1.0, 0.0, 1.2), 1, pytest.approx(0.2, abs=1e-9)),
        ((1.0, 0.0, 1.2), 2, 0),
        ((1.0, 1.0, 1.5), 1, 0),
        ((0.0, 1.0, 1.4), 1, pytest.approx(SQRT2 - 1.4, abs=1e-9)),
        ((-1.0, -1.0, -1.0), 1, 0),
    ],
)
def test_eval2_specific_values(square_landscape, grade, k, expected):
    """Re-runs the eval-spec table through eval2 on a single-row grade array."""
    out = square_landscape.eval2(np.array([grade]), k)
    assert out[0, k - 1] == expected


def test_eval2_empty_landscape_returns_zeros(triangle_landscape):
    grades = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [-1.0, 5.0, 0.5]])
    out = triangle_landscape.eval2(grades, 4)
    np.testing.assert_array_equal(out, np.zeros((3, 4)))


def test_eval2_grades_below_all_births_are_zero(square_landscape):
    """No birth ≤ grade → no contributions → all zeros."""
    grades = np.array([[-1.0, -1.0, -1.0], [0.5, -0.5, 0.0]])
    out = square_landscape.eval2(grades, 3)
    np.testing.assert_array_equal(out, np.zeros((2, 3)))


def test_eval2_above_layer_count_is_zero(square_landscape):
    """Square has at most 2 simultaneous bars; layers 3+ must be 0."""
    grades = np.array([[1.0, 1.0, 5.0], [1.0, 0.0, 1.2]])
    out = square_landscape.eval2(grades, 5)
    np.testing.assert_array_equal(out[:, 2:], np.zeros((2, 3)))


def test_eval2_layers_monotone_decreasing(random_landscape):
    """Layer k+1 ≤ layer k everywhere by definition of a persistence landscape."""
    grades = _grid_around_landscape(random_landscape, n_per_axis=5)
    out = random_landscape.eval2(grades, 4)
    diffs = np.diff(out, axis=1)
    assert (diffs <= 1e-12).all(), f"non-monotone layers found: max excess = {diffs.max()}"


def test_eval2_deterministic(square_landscape):
    grades = _grid_around_landscape(square_landscape, n_per_axis=4)
    a = square_landscape.eval2(grades, 3)
    b = square_landscape.eval2(grades, 3)
    np.testing.assert_array_equal(a, b)


def test_eval2_grade_order_independent(random_landscape):
    """Permuting the input grades just permutes the rows of the output."""
    grades = _grid_around_landscape(random_landscape, n_per_axis=4)
    perm = np.random.default_rng(0).permutation(len(grades))
    base = random_landscape.eval2(grades, 3)
    permuted = random_landscape.eval2(grades[perm], 3)
    np.testing.assert_array_equal(permuted, base[perm])


def test_eval2_single_cell(square_landscape):
    out = square_landscape.eval2(np.array([[1.0, 0.0, 1.2]]), 2)
    assert out.shape == (1, 2)
    assert out[0, 0] == pytest.approx(0.2, abs=1e-9)
    assert out[0, 1] == 0.0


def test_eval2_k_max_zero_returns_empty_columns(square_landscape):
    grades = np.array([[1.0, 0.0, 1.2]])
    out = square_landscape.eval2(grades, 0)
    assert out.shape == (1, 0)


def test_eval2_rejects_negative_k_max(square_landscape):
    with pytest.raises(ValueError):
        square_landscape.eval2(np.array([[1.0, 0.0, 1.2]]), -1)


def test_eval2_rejects_non_2d_grades(square_landscape):
    with pytest.raises(ValueError):
        square_landscape.eval2(np.array([1.0, 0.0, 1.2]), 1)


@pytest.mark.parametrize("fixture_name", list(COMPLEX_FIXTURES))
def test_eval2_matches_eval_on_complex_fixtures(fixture_name):
    """Same cross-check on the complex fixtures. Grid is intentionally tiny so the
    O(cells * pairings * syzygies) Python reference doesn't dominate the suite."""
    trajectories, max_metric_value, hom_dim = COMPLEX_FIXTURES[fixture_name]()
    landscape = compute_spatiotemporal_landscapes_sparse(trajectories, max_metric_value, hom_dim)
    grades = _grid_around_landscape(landscape, n_per_axis=4)
    k_max = 3
    out = landscape.eval2(grades, k_max)
    np.testing.assert_allclose(out, _eval_cellwise(landscape, grades, k_max), atol=1e-12)


# ---------- tests: save / load roundtrip ----------

def test_save_load_roundtrip(square_landscape):
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "landscape.slan")
        square_landscape.save(path, precision=12)
        loaded = SparseLandscape.load(path)
    _assert_pairings_equal(loaded.pairings, square_landscape.pairings, tol=1e-10)


# ---------- tests: Grade helpers (used by SparseLandscape.eval) ----------

class TestGrade:

    def test_le_componentwise(self):
        assert Grade((1, 2, 3)) <= Grade((1, 2, 3))
        assert Grade((1, 2, 3)) <= Grade((1, 3, 3))
        assert not (Grade((1, 2, 3)) <= Grade((0, 2, 3)))

    def test_lt_strict(self):
        # equal is not strictly less
        assert not (Grade((1, 2, 3)) < Grade((1, 2, 3)))
        # ≤ in all + strictly less in at least one
        assert Grade((1, 2, 3)) < Grade((1, 2, 4))
        # incomparable
        assert not (Grade((1, 2, 3)) < Grade((0, 3, 3)))

    def test_ge_gt_dual(self):
        assert Grade((1, 2, 4)) > Grade((1, 2, 3))
        assert not (Grade((1, 2, 3)) > Grade((1, 2, 3)))
        assert Grade((1, 2, 3)) >= Grade((1, 2, 3))

    def test_compare_size_mismatch_raises(self):
        with pytest.raises(ValueError):
            _ = Grade((1, 2)) <= Grade((1, 2, 3))

    def test_join_takes_componentwise_max(self):
        assert Grade.join(Grade((1, 5, 2)), Grade((4, 3, 2))) == Grade((4, 5, 2))

    def test_colex_compare_orders_by_last_first(self):
        # last component dominates colex order
        assert Grade.colex_compare(Grade((9, 9, 1)), Grade((0, 0, 2))) == -1
        # tie on last → next-to-last decides
        assert Grade.colex_compare(Grade((9, 1, 5)), Grade((0, 2, 5))) == -1
        # equal
        assert Grade.colex_compare(Grade((1, 2, 3)), Grade((1, 2, 3))) == 0

    def test_slice_returns_grade(self):
        g = Grade((1, 2, 3))
        assert isinstance(g[:-1], Grade)
        assert g[:-1] == Grade((1, 2))


# ---------- complex (~1s) fixtures: snapshot regression ----------
#
# Snapshots live at test/snapshots/<name>.json and are regenerated by running
# `python test/update_complex_snapshots.py`. Don't edit them by hand — failures
# here mean either a behavioural regression or a deliberate change that
# warrants a snapshot refresh.

@pytest.mark.parametrize("fixture_name", list(COMPLEX_FIXTURES))
def test_complex_fixture_snapshot(fixture_name):
    snapshot_path = SNAPSHOT_DIR / f"{fixture_name}.json"
    assert snapshot_path.exists(), (
        f"missing snapshot {snapshot_path}. Run "
        "`python test/update_complex_snapshots.py` to generate it."
    )

    snapshot = json.loads(snapshot_path.read_text())
    trajectories, max_metric_value, hom_dim = COMPLEX_FIXTURES[fixture_name]()

    # Sanity-check the fixture matches the snapshot's recorded inputs — guards
    # against the case where someone edits a fixture but forgets to regenerate.
    assert list(trajectories.shape) == snapshot["trajectory_shape"]
    assert max_metric_value == snapshot["max_metric_value"]
    assert hom_dim == snapshot["hom_dim"]

    res = compute_spatiotemporal_landscapes_sparse(trajectories, max_metric_value, hom_dim)
    expected = [(tuple(b), [tuple(s) for s in syz]) for b, syz in snapshot["pairings"]]
    _assert_pairings_equal(res.pairings, expected)
