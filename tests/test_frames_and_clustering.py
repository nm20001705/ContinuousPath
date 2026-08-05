"""Tests for the (u,v) basis and the edge-clustering that finds solid pieces.

basis_vectors defines the frame that BOTH the cached wing slices and the
bridge/hole placement are expressed in. If those two ever disagree the
features land in the wrong plane, with no error anywhere -- which is why
the definition was consolidated into one function. These tests pin the
properties that consolidation relies on.
"""

import numpy as np
import pytest

from rib_slice_core import basis_vectors, cluster_segments, prepare_segment


AXES = {
    'X': np.array([1.0, 0.0, 0.0]),
    'Y': np.array([0.0, 1.0, 0.0]),
    'Z': np.array([0.0, 0.0, 1.0]),
}


# ---------------------------------------------------------------------
# basis_vectors
# ---------------------------------------------------------------------

@pytest.mark.parametrize('name', sorted(AXES))
def test_basis_is_orthonormal_for_every_principal_axis(name):
    """The `abs(prim[0]) > 0.9` branch exists to dodge a degenerate cross
    product when prim is X. Both branches must produce a proper frame."""
    prim, u, v = basis_vectors(AXES[name])

    for vec, label in ((prim, 'prim'), (u, 'u'), (v, 'v')):
        assert np.linalg.norm(vec) == pytest.approx(1.0), f"{label} not unit"

    assert np.dot(prim, u) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(prim, v) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(u, v) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize('name', sorted(AXES))
def test_basis_spans_the_plane_perpendicular_to_prim(name):
    """u x v must recover prim (up to sign) -- i.e. the pair really does
    span the slicing plane rather than some skew subspace."""
    prim, u, v = basis_vectors(AXES[name])
    recovered = np.cross(u, v)
    assert abs(abs(np.dot(recovered, prim)) - 1.0) < 1e-12


def test_basis_normalises_an_unnormalised_direction():
    """Callers pass primary_dir straight from config, which need not be
    unit length."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 7.5]))
    assert np.linalg.norm(prim) == pytest.approx(1.0)
    assert prim[2] == pytest.approx(1.0)


def test_basis_is_deterministic():
    """Repeated calls must agree exactly -- the cached slices on disk are
    only reusable if the frame is reproducible across runs."""
    a = basis_vectors(AXES['Z'])
    b = basis_vectors(AXES['Z'])
    for x, y in zip(a, b):
        assert np.array_equal(x, y)


def test_basis_handles_an_arbitrary_oblique_direction():
    prim, u, v = basis_vectors(np.array([0.3, -0.5, 0.81]))
    assert np.linalg.norm(prim) == pytest.approx(1.0)
    assert np.dot(prim, u) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(u, v) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------
# cluster_segments
# ---------------------------------------------------------------------

def test_two_disjoint_pieces_stay_separate():
    """A slice through a rib near a cutout yields several disconnected
    chords; keeping them apart is what lets iter_solid_pieces pick the one
    the rib centreline actually passes through."""
    a = np.array([[0.0, 0.0], [1.0, 0.0], [10.0, 0.0], [11.0, 0.0]])
    b = np.array([[1.0, 0.0], [2.0, 0.0], [11.0, 0.0], [12.0, 0.0]])
    groups = cluster_segments(a, b)
    assert len(groups) == 2
    assert sorted(sorted(g) for g in groups) == [[0, 1], [2, 3]]


def test_segments_sharing_an_endpoint_merge_into_one_piece():
    a = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    b = np.array([[1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    groups = cluster_segments(a, b)
    assert len(groups) == 1
    assert sorted(groups[0]) == [0, 1, 2]


def test_a_closed_loop_is_a_single_piece():
    """A closed chord ring must not be split at its seam."""
    a = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    b = np.array([[1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]])
    groups = cluster_segments(a, b)
    assert len(groups) == 1
    assert sorted(groups[0]) == [0, 1, 2, 3]


def test_every_segment_is_assigned_exactly_once():
    """No segment may be dropped or double-counted -- a lost segment
    silently shrinks the measured solid extent."""
    rng = np.random.default_rng(0)
    a = rng.random((20, 2))
    b = rng.random((20, 2))
    groups = cluster_segments(a, b)
    flat = sorted(i for g in groups for i in g)
    assert flat == list(range(20))


def test_clustering_joins_endpoints_within_rounding_tolerance():
    """Endpoints are matched on rounded coordinates, so vertices that
    agree to the tolerance are treated as shared."""
    a = np.array([[0.0, 0.0], [1.0, 0.0]])
    b = np.array([[1.0, 0.0], [1.0 + 1e-9, 0.0]])
    groups = cluster_segments(a, b, tol=5)
    assert len(groups) == 1


def test_single_segment_is_one_group():
    groups = cluster_segments(np.array([[0.0, 0.0]]), np.array([[1.0, 0.0]]))
    assert groups == [[0]]


# ---------------------------------------------------------------------
# prepare_segment
# ---------------------------------------------------------------------

def test_prepare_segment_builds_an_in_plane_frame():
    """line_dir_3d must be perpendicular to both the slab normal and the
    slicing direction -- it is the axis bridge/hole widths are measured
    along."""
    prim, u, v = basis_vectors(AXES['Z'])
    seg = {
        'p0': np.array([0.0, 0.0, 0.0]),
        'dir': np.array([1.0, 0.0, 0.0]),
        'slab_normal': np.array([0.0, 1.0, 0.0]),
    }
    prep = prepare_segment(seg, prim, u, v)
    assert prep is not None
    assert np.dot(prep['line_dir_3d'], seg['slab_normal']) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(prep['line_dir_3d'], prim) == pytest.approx(0.0, abs=1e-12)
    assert np.linalg.norm(prep['line_dir_3d']) == pytest.approx(1.0)
    assert np.linalg.norm(prep['bridge_dir_2d']) == pytest.approx(1.0)


def test_prepare_segment_rejects_a_rib_plane_perpendicular_to_prim():
    """When the slab normal is parallel to the slicing direction the rib
    plane is perpendicular to it, cross() collapses, and no meaningful
    width axis exists. Returning None is what makes that rib get skipped
    rather than produce garbage -- this is the degenerate case that occurs
    when (grid_orientation +/- rib_angle) hits 90 degrees."""
    prim, u, v = basis_vectors(AXES['Z'])
    seg = {
        'p0': np.array([0.0, 0.0, 0.0]),
        'dir': np.array([1.0, 0.0, 0.0]),
        'slab_normal': np.array([0.0, 0.0, 1.0]),   # parallel to prim
    }
    assert prepare_segment(seg, prim, u, v) is None


def test_prepare_segment_plane_offset_is_invariant_along_the_rib():
    """p0 is only ever 'some point on the rib line'. Moving it along the
    rib direction must not change the derived plane -- several callers
    rely on this, and it is why a mismatched segment_bounds entry was
    harmless."""
    prim, u, v = basis_vectors(AXES['Z'])
    base = {
        'dir': np.array([1.0, 0.0, 0.0]),
        'slab_normal': np.array([0.0, 1.0, 0.0]),
    }
    a = prepare_segment({**base, 'p0': np.array([0.0, 3.0, 0.0])}, prim, u, v)
    b = prepare_segment({**base, 'p0': np.array([50.0, 3.0, 0.0])}, prim, u, v)
    assert a['plane_offset'] == pytest.approx(b['plane_offset'])
