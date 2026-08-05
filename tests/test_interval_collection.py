"""Tests for collect_line_intervals, which turns per-slice width decisions
into the (d0, d1, t_lo, t_hi) strips that get extruded into solids.

The subtle part is that t values are stored RELATIVE to t_rib -- the point
where the rib centreline crosses that height. Storing raw t instead adds a
height-dependent shift to every strip, which showed up as bridges and
holes drifting sideways as they climbed. Pinned here because the symptom
(drift) is far from the cause (one missing subtraction).
"""

import numpy as np
import pytest
import trimesh

from rib_slice_core import basis_vectors, collect_line_intervals


def slab_mesh(length=20.0, height=10.0, thickness=0.2):
    """A flat rib-segment slab lying in the x-z plane, centred on x=0."""
    m = trimesh.creation.box(extents=(length, thickness, height))
    m.apply_translation((0.0, 0.0, height / 2.0))
    return m


def straight_seg():
    """A rib line along x, with its slab plane normal along y."""
    return {
        'p0': np.array([0.0, 0.0, 0.0]),
        'dir': np.array([0.0, 0.0, 1.0]),
        'slab_normal': np.array([0.0, 1.0, 0.0]),
    }


def full_width(piece):
    return piece['t_min_solid'], piece['t_max_solid']


def test_intervals_pair_consecutive_slices():
    """N usable slices produce N-1 strips, each spanning one z step."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    z_vals = np.arange(1.0, 9.0, 1.0)

    intervals, frame = collect_line_intervals(
        slab_mesh(), straight_seg(), prim, u, v, z_vals, full_width)

    assert frame is not None
    assert len(intervals) >= 1
    for d0, d1, t_lo, t_hi in intervals:
        assert d1 > d0
        assert t_hi >= t_lo
    ds = [iv[0] for iv in intervals]
    assert ds == sorted(ds), "intervals must be ordered by height"


def test_intervals_are_empty_when_policy_rejects_everything():
    """A policy that never returns a width yields no strips at all,
    rather than a degenerate one."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    z_vals = np.arange(1.0, 9.0, 1.0)

    intervals, _ = collect_line_intervals(
        slab_mesh(), straight_seg(), prim, u, v, z_vals, lambda p: None)

    assert intervals == []


def test_fewer_than_two_slices_yields_no_intervals():
    """One slice cannot form a strip -- a strip needs a start and an end."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    calls = {'n': 0}

    def only_first(piece):
        calls['n'] += 1
        return full_width(piece) if calls['n'] == 1 else None

    intervals, frame = collect_line_intervals(
        slab_mesh(), straight_seg(), prim, u, v,
        np.arange(1.0, 9.0, 1.0), only_first)

    assert intervals == []
    assert frame is not None, "the frame is still reported for the one slice"


def test_t_values_are_stored_relative_to_the_centreline():
    """The regression guard: a slab centred on the rib centreline must
    produce strips centred on ZERO, because t is stored as an offset from
    t_rib. If raw t were stored, the offsets would inherit the slab's
    absolute position and drift with height."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    z_vals = np.arange(1.0, 9.0, 1.0)

    intervals, _ = collect_line_intervals(
        slab_mesh(length=20.0), straight_seg(), prim, u, v, z_vals, full_width)

    assert intervals
    for _, _, t_lo, t_hi in intervals:
        centre = (t_lo + t_hi) / 2.0
        assert centre == pytest.approx(0.0, abs=1e-6), (
            "strips must be centred on the rib centreline, not on the "
            "global origin")


def test_offsets_do_not_drift_with_height():
    """Every slice of a prismatic slab has the same cross-section, so the
    stored offsets must be identical at every height. Drift here is
    exactly the bug this test exists for."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    z_vals = np.arange(1.0, 9.0, 1.0)

    intervals, _ = collect_line_intervals(
        slab_mesh(), straight_seg(), prim, u, v, z_vals, full_width)

    los = [round(iv[2], 6) for iv in intervals]
    his = [round(iv[3], 6) for iv in intervals]
    assert len(set(los)) == 1, f"t_lo drifts with height: {sorted(set(los))}"
    assert len(set(his)) == 1, f"t_hi drifts with height: {sorted(set(his))}"


def tilted_seg(tilt=0.3):
    """A rib line leaning in x as it climbs in z.

    This is the case that matters: with an upright rib the centreline
    crossing sits at t_rib = 0 for every slice, so storing raw t and
    storing t - t_rib are indistinguishable. Only a tilt makes t_rib
    non-zero AND height-dependent, which is exactly the situation the
    drift bug corrupted.
    """
    d = np.array([tilt, 0.0, 1.0])
    return {
        'p0': np.array([0.0, 0.0, 0.0]),
        'dir': d / np.linalg.norm(d),
        'slab_normal': np.array([0.0, 1.0, 0.0]),
    }


def test_stored_offsets_reconstruct_onto_the_real_geometry():
    """The strongest form of the drift guard.

    solidify_rib_line rebuilds each point as
        origin_const + d*axis_d + t*line_dir_3d
    so a correctly stored offset must land back on the slab it came from.
    Storing raw t instead leaves in the t_rib term, displacing every strip
    sideways by an amount that grows with height -- points then land
    outside the slab entirely. Checking the reconstruction, rather than
    the stored numbers, catches that regardless of how the offset is
    represented internally.
    """
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    z_vals = np.arange(1.0, 9.0, 1.0)
    half_length = 10.0
    mesh = slab_mesh(length=2 * half_length)

    intervals, frame = collect_line_intervals(
        mesh, tilted_seg(), prim, u, v, z_vals, full_width)

    assert intervals, "expected strips"
    origin_const = frame['origin_const']
    axis_d = frame['axis_d']
    line_dir = frame['line_dir_3d']
    tilt = 0.3

    for d0, d1, t_lo, t_hi in intervals:
        # Each strip spans a pair of slices and stores the bounding range
        # of their four t-values (the documented conservative choice over
        # an exact trapezoid). On a tilted rib the centreline moves by
        # tilt*(d1-d0) between them, so reconstructing at d0 may legally
        # overshoot the slab by exactly that much -- and no more. The
        # drift bug's error instead grows with absolute height, far
        # exceeding one step.
        allowed = half_length + tilt * (d1 - d0) + 1e-6
        for t in (t_lo, t_hi):
            pt = origin_const + d0 * axis_d + t * line_dir
            assert abs(pt[0]) <= allowed, (
                f"reconstructed x={pt[0]:.3f} at height {d0} exceeds the "
                f"slab half-length {half_length} by more than one step's "
                f"drift -- offsets are not relative to the rib centreline")
            assert pt[2] == pytest.approx(d0, abs=1e-6), (
                "reconstructed height must match the slice it came from")


def test_offsets_track_the_centreline_on_a_tilted_rib():
    """A prismatic slab has the same absolute t-extent at every height, so
    once the moving t_rib is subtracted the stored offsets MUST shift with
    height. Constant offsets here mean t_rib was never removed."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    z_vals = np.arange(1.0, 9.0, 1.0)

    intervals, _ = collect_line_intervals(
        slab_mesh(length=20.0), tilted_seg(), prim, u, v, z_vals, full_width)

    los = [round(iv[2], 6) for iv in intervals]
    assert len(set(los)) > 1, (
        "offsets are constant on a tilted rib -- t_rib was not subtracted")


def test_frame_carries_the_axes_needed_to_rebuild_3d():
    """solidify_rib_line reconstructs points as
    origin_const + d*axis_d + t*line_dir_3d, so all three must be present
    and the direction vectors must be unit length."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    _, frame = collect_line_intervals(
        slab_mesh(), straight_seg(), prim, u, v,
        np.arange(1.0, 9.0, 1.0), full_width)

    assert frame is not None
    for key in ('origin_const', 'axis_d', 'line_dir_3d', 'slab_normal'):
        assert key in frame, f"frame missing {key}"
    assert np.linalg.norm(frame['line_dir_3d']) == pytest.approx(1.0)
    assert np.linalg.norm(frame['slab_normal']) == pytest.approx(1.0)


def test_strip_bounds_enclose_both_endpoints_of_each_pair():
    """Each strip uses the bounding range of its two slices' t-values --
    a deliberately conservative choice over the exact trapezoid."""
    prim, u, v = basis_vectors(np.array([0.0, 0.0, 1.0]))
    widths = iter([(-1.0, 1.0), (-3.0, 0.5), (-2.0, 2.0)])
    seen = []

    def varying(piece):
        try:
            w = next(widths)
        except StopIteration:
            return None
        seen.append((w[0] - piece['t_rib'], w[1] - piece['t_rib']))
        return w

    intervals, _ = collect_line_intervals(
        slab_mesh(), straight_seg(), prim, u, v,
        np.arange(1.0, 9.0, 1.0), varying)

    assert len(intervals) == len(seen) - 1
    for k, (_, _, t_lo, t_hi) in enumerate(intervals):
        a, b = seen[k], seen[k + 1]
        assert t_lo == pytest.approx(min(a[0], a[1], b[0], b[1]))
        assert t_hi == pytest.approx(max(a[0], a[1], b[0], b[1]))
