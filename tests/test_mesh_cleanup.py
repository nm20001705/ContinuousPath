"""Tests for post-boolean mesh cleanup.

drop_sliver_components removes the degenerate scraps a boolean sheds
along 0.1mm slit edges. Getting this wrong is unusually easy: the first
two implementations each destroyed the mesh in a different way, and both
looked fine until the exported STL was re-examined. The guard those
failures produced is pinned here.
"""

import numpy as np
import pytest
import trimesh

from mesh_simplify_utils import drop_sliver_components


def box(extents=(1.0, 1.0, 1.0), translate=(0.0, 0.0, 0.0)):
    m = trimesh.creation.box(extents=extents)
    m.apply_translation(translate)
    return m


def zero_volume_shell(translate=(10.0, 0.0, 0.0)):
    """A closed but perfectly flat shell -- two coincident triangles.
    This is what the boolean actually leaves behind: watertight in the
    edge-count sense, with |volume| == 0."""
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    f = np.array([[0, 1, 2], [0, 2, 1]])
    m = trimesh.Trimesh(vertices=v, faces=f, process=False)
    m.apply_translation(translate)
    return m


def test_zero_volume_scraps_are_removed():
    mesh = trimesh.util.concatenate([box(), zero_volume_shell()])
    assert len(mesh.split(only_watertight=False)) == 2

    out = drop_sliver_components(mesh, verbose=False)

    assert len(out.split(only_watertight=False)) == 1
    assert out.volume == pytest.approx(1.0, rel=1e-6)


def test_solid_volume_is_preserved_exactly():
    """Removing zero-volume debris must not change the part's volume --
    this is the check that proved the real cleanup was lossless."""
    mesh = trimesh.util.concatenate(
        [box(), zero_volume_shell((10, 0, 0)), zero_volume_shell((20, 0, 0))])
    before = mesh.volume
    out = drop_sliver_components(mesh, verbose=False)
    assert out.volume == pytest.approx(before, rel=1e-12)


def test_genuinely_detached_solid_pieces_are_kept():
    """Full-depth rib slabs can legitimately isolate a cell, so a second
    real solid must survive. Dropping by *volume* rather than by size or
    component count is what makes that safe."""
    mesh = trimesh.util.concatenate([
        box(extents=(2.0, 2.0, 2.0)),
        box(extents=(0.5, 0.5, 0.5), translate=(10.0, 0.0, 0.0)),
        zero_volume_shell((20.0, 0.0, 0.0)),
    ])
    out = drop_sliver_components(mesh, verbose=False)
    comps = out.split(only_watertight=False)
    assert len(comps) == 2, "the small but solid piece must be kept"
    assert out.volume == pytest.approx(8.0 + 0.125, rel=1e-6)


def test_small_but_solid_piece_survives_when_above_threshold():
    """min_volume is a degeneracy floor, not a size filter."""
    mesh = trimesh.util.concatenate([
        box(extents=(1.0, 1.0, 1.0)),
        box(extents=(0.1, 0.1, 0.1), translate=(5.0, 0.0, 0.0)),  # vol 1e-3
    ])
    out = drop_sliver_components(mesh, min_volume=1e-6, verbose=False)
    assert len(out.split(only_watertight=False)) == 2


def flat_ribbon(translate=(10.0, 0.0, 0.0)):
    """A long, flat, effectively zero-volume ribbon with MANY faces.

    Modelled on the real debris (e.g. 10.55 x 0.09 x 0.00 mm). Built as a
    subdivided near-degenerate box so it is a single connected component,
    unlike a double-sided sheet whose 4-face edges break adjacency.

    The point is that it has MORE faces than a small *solid* piece, so a
    cleanup filtering on face count or bounding size would keep this junk
    while discarding real material. Only volume separates them.
    """
    m = trimesh.creation.box(extents=(10.0, 0.09, 1e-9))
    m = m.subdivide()
    m.apply_translation(translate)
    return m


def test_many_faced_zero_volume_debris_is_still_dropped():
    """Guards the choice of criterion: this ribbon has far more faces than
    the small solid cube below, so filtering on face count or bounding
    size would keep the debris and discard real material. Only volume
    separates them correctly."""
    ribbon = flat_ribbon()
    small_solid = box(extents=(0.4, 0.4, 0.4), translate=(30.0, 0.0, 0.0))
    assert len(ribbon.faces) > len(small_solid.faces)

    mesh = trimesh.util.concatenate([box(extents=(2.0, 2.0, 2.0)),
                                     ribbon, small_solid])
    out = drop_sliver_components(mesh, verbose=False)
    comps = out.split(only_watertight=False)

    assert len(comps) == 2, "ribbon dropped, both solids kept"
    assert out.volume == pytest.approx(8.0 + 0.064, rel=1e-6)


def test_single_component_mesh_is_returned_untouched():
    mesh = box()
    out = drop_sliver_components(mesh, verbose=False)
    assert out is mesh


def test_none_passes_through():
    assert drop_sliver_components(None, verbose=False) is None


def test_mesh_with_no_slivers_is_returned_unchanged():
    """If nothing would be dropped, hand back the original object rather
    than a rebuilt copy."""
    mesh = trimesh.util.concatenate([box(), box(translate=(5.0, 0.0, 0.0))])
    out = drop_sliver_components(mesh, verbose=False)
    assert out is mesh


def test_watertightness_is_never_traded_away():
    """The guard that caught the first broken implementation: if removing
    components would leave the kept geometry open, keep the mesh intact
    instead. A leaky cut tool sews into a shell, Part.makeSolid() then
    fails, and the wing cut returns a Null shape or runs for hours -- so
    silently returning something leaky is far worse than keeping debris."""
    mesh = trimesh.util.concatenate([box(), zero_volume_shell()])
    assert mesh.is_watertight

    out = drop_sliver_components(mesh, verbose=False)

    assert out.is_watertight, "cleanup must not open a watertight mesh"


def test_everything_below_threshold_keeps_the_mesh_intact():
    """If no component qualifies, returning an empty mesh would silently
    delete the part."""
    mesh = trimesh.util.concatenate(
        [zero_volume_shell((0, 0, 0)), zero_volume_shell((10, 0, 0))])
    out = drop_sliver_components(mesh, verbose=False)
    assert len(out.faces) == len(mesh.faces)
