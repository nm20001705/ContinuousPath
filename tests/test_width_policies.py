"""Tests for the bridge and hole width policies in rib_slice_core.

These two functions decide, for every Z-slice of every rib, how wide a
strip of material to cut. They are pure arithmetic on a `piece` dict, so
they need neither FreeCAD nor real geometry -- yet both regressed
silently during development and each cost hours of full pipeline runs to
diagnose. That asymmetry is the whole reason these tests exist.
"""

import numpy as np
import pytest

from rib_slice_core import bridge_width_interval, hole_width_interval


def piece(t_min, t_max, d=5.0, d_min=0.0, d_max=10.0):
    """Minimal `piece` dict -- only the keys the width policies read."""
    return {
        't_min_solid': t_min,
        't_max_solid': t_max,
        'd': d,
        'd_min': d_min,
        'd_max': d_max,
    }


# ---------------------------------------------------------------------
# bridge_width_interval
# ---------------------------------------------------------------------

def test_bridge_uses_full_height_when_it_fits():
    """A wide section gets exactly bridge_height, centred."""
    lo, hi = bridge_width_interval(piece(0.0, 10.0), bridge_height=1.0)
    assert hi - lo == pytest.approx(1.0)
    assert (lo + hi) / 2 == pytest.approx(5.0)


def test_bridge_is_centred_on_the_piece_not_the_origin():
    """The bridge centres on the solid extent, wherever that sits."""
    lo, hi = bridge_width_interval(piece(100.0, 110.0), bridge_height=2.0)
    assert (lo + hi) / 2 == pytest.approx(105.0)
    assert hi - lo == pytest.approx(2.0)


def test_bridge_clamps_instead_of_skipping_when_section_is_thin():
    """The regression this session: a section narrower than bridge_height
    used to return None, which broke the bridge column into disconnected
    pieces exactly where the wall is thinnest -- near the leading/trailing
    edge and the tip, i.e. where the bridge matters most. It must now
    clamp to the available width instead."""
    lo, hi = bridge_width_interval(piece(0.0, 0.35), bridge_height=1.0)
    assert (lo, hi) != (None, None)
    assert hi - lo == pytest.approx(0.35), "should fill the section, not skip it"
    assert lo == pytest.approx(0.0)
    assert hi == pytest.approx(0.35)


def test_bridge_never_exceeds_the_available_solid():
    """Clamping must not spill outside the piece -- extruding beyond the
    solid would put bridge material where there is no rib to attach to."""
    for width in (0.05, 0.5, 0.99, 1.0, 5.0):
        lo, hi = bridge_width_interval(piece(2.0, 3.0), bridge_height=width)
        assert lo >= 2.0 - 1e-12
        assert hi <= 3.0 + 1e-12


def test_bridge_returns_none_only_on_degenerate_extent():
    """A zero-width strip is not a thin bridge -- it is geometry that
    would break the extrude and the downstream boolean."""
    assert bridge_width_interval(piece(5.0, 5.0), bridge_height=1.0) is None
    assert bridge_width_interval(piece(5.0, 4.0), bridge_height=1.0) is None
    assert bridge_width_interval(piece(0.0, 1e-9), bridge_height=1.0) is None


def test_bridge_accepts_a_section_at_the_min_width_threshold():
    """Just above the degeneracy floor is still a real bridge."""
    result = bridge_width_interval(piece(0.0, 1e-5), bridge_height=1.0,
                                   min_width=1e-6)
    assert result is not None
    lo, hi = result
    assert hi - lo == pytest.approx(1e-5)


# ---------------------------------------------------------------------
# hole_width_interval
# ---------------------------------------------------------------------

def always_full(x):
    return 1.0


def test_hole_margin_includes_half_the_wall_thickness():
    """eff_margin is hole_margin + thickness/2, not hole_margin alone --
    the hole must clear the neighbouring rib's *extruded* volume, not just
    its centre plane."""
    p = piece(0.0, 10.0, d=5.0, d_min=0.0, d_max=10.0)
    lo, hi = hole_width_interval(p, always_full, hole_margin=1.0, thickness=2.0)
    # eff_margin = 1.0 + 1.0 = 2.0, so the usable span is [2, 8]
    assert lo == pytest.approx(2.0)
    assert hi == pytest.approx(8.0)


def test_hole_condition_is_mapped_onto_the_margin_free_span():
    """x=0 must land at d_min+eff_margin and x=1 at d_max-eff_margin, NOT
    at the raw geometric bounds. A profile that tapers to zero at its ends
    should therefore vanish exactly at the inset boundary."""
    seen = []

    def record(x):
        seen.append(x)
        return 1.0

    # eff_margin = 0.5, so the usable d-span is [0.5, 9.5]
    hole_width_interval(piece(0.0, 10.0, d=0.5), record,
                        hole_margin=0.5, thickness=0.0)
    hole_width_interval(piece(0.0, 10.0, d=9.5), record,
                        hole_margin=0.5, thickness=0.0)
    hole_width_interval(piece(0.0, 10.0, d=5.0), record,
                        hole_margin=0.5, thickness=0.0)

    assert seen[0] == pytest.approx(0.0), "x=0 should be the inset start"
    assert seen[1] == pytest.approx(1.0), "x=1 should be the inset end"
    assert seen[2] == pytest.approx(0.5), "midpoint should be x=0.5"


def test_hole_absent_inside_the_d_margin_band():
    """Slices within eff_margin of the top or bottom get no hole at all --
    the point_condition is never consulted there."""
    called = []

    def record(x):
        called.append(x)
        return 1.0

    for d in (0.0, 0.4, 9.6, 10.0):
        p = piece(0.0, 10.0, d=d, d_min=0.0, d_max=10.0)
        assert hole_width_interval(p, record, hole_margin=0.5,
                                   thickness=0.0) is None
    assert called == [], "point_condition must not be consulted in the margin"


def test_hole_none_when_section_too_thin_for_the_margin():
    """If the two t-margins overlap there is no room for a hole."""
    p = piece(0.0, 1.0, d=5.0, d_min=0.0, d_max=10.0)
    assert hole_width_interval(p, always_full, hole_margin=1.0,
                               thickness=0.0) is None


def test_hole_width_scales_with_the_condition_value():
    """A factor of f yields f times the available width, centred."""
    p = piece(0.0, 10.0, d=5.0, d_min=0.0, d_max=10.0)
    lo, hi = hole_width_interval(p, lambda x: 0.5, hole_margin=0.0,
                                 thickness=0.0)
    assert hi - lo == pytest.approx(5.0)
    assert (lo + hi) / 2 == pytest.approx(5.0)


def test_hole_never_exceeds_the_margin_free_span():
    """Even a condition greater than 1 must stay inside the clearance."""
    p = piece(0.0, 10.0, d=5.0, d_min=0.0, d_max=10.0)
    lo, hi = hole_width_interval(p, lambda x: 5.0, hole_margin=1.0,
                                 thickness=0.0)
    assert lo >= 1.0 - 1e-12
    assert hi <= 9.0 + 1e-12


def test_hole_none_when_condition_collapses_to_zero():
    """A zero-width hole is degenerate geometry, not a tiny hole."""
    p = piece(0.0, 10.0, d=5.0, d_min=0.0, d_max=10.0)
    assert hole_width_interval(p, lambda x: 0.0, hole_margin=0.0,
                               thickness=0.0) is None


def test_circular_condition_is_blunt_at_the_ends_not_pointed():
    """The semicircle profile has a vertical tangent at x=0 and x=1, so
    the hole opens abruptly rather than tapering to a point. This looked
    like a bug during development but is inherent to the shape -- pinned
    here so the observation is not re-litigated."""
    def circle(x):
        return np.sqrt(max(0.0, 1 - (2 * x - 1) ** 2))

    p_near_end = piece(0.0, 10.0, d=0.05, d_min=0.0, d_max=10.0)
    p_middle = piece(0.0, 10.0, d=5.0, d_min=0.0, d_max=10.0)

    near = hole_width_interval(p_near_end, circle, 0.0, 0.0)
    mid = hole_width_interval(p_middle, circle, 0.0, 0.0)

    near_w = near[1] - near[0]
    mid_w = mid[1] - mid[0]
    # 0.5% along the span already reaches >13% of full width
    assert near_w / mid_w > 0.13
