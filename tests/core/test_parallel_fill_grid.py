"""
Tests for optimized fill_grid implementation
"""

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy import units as u

from hifast.core.parallel_fill_grid import optimized_search_around_sky


def test_optimized_vs_original():
    """Test that optimized version produces same results as original"""
    np.random.seed(42)
    n_cata = 100_000
    n_grid = 20_000
    
    cata = SkyCoord(
        np.random.uniform(0, 30, n_cata),
        np.random.uniform(20, 40, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(10, 20, n_grid),
        np.random.uniform(25, 35, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 90 * u.arcsec
    
    # Run both versions
    # Note: cata.search_around_sky(grid) internally calls search_around_sky(grid, cata)
    # So we need to call optimized_search_around_sky(grid, cata) to match
    idx1, idx2, _, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=4, verbose=False
    )
    idx1_orig, idx2_orig, _, _ = cata.search_around_sky(grid, r_cut)
    
    # Compare results
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"Results mismatch: {len(pairs_opt ^ pairs_orig)} pairs differ"
    assert len(idx1) == len(idx1_orig), "Different number of matches"


def test_optimized_single_worker():
    """Test optimized version with single worker"""
    np.random.seed(123)
    n_cata = 10_000
    n_grid = 5_000
    
    cata = SkyCoord(
        np.random.uniform(0, 10, n_cata),
        np.random.uniform(0, 10, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(2, 8, n_grid),
        np.random.uniform(2, 8, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 60 * u.arcsec
    
    idx1, idx2, d2d, d3d = optimized_search_around_sky(
        cata, grid, r_cut, n_workers=1, verbose=False
    )
    
    assert len(idx1) > 0, "Should find some matches"
    assert len(idx1) == len(idx2), "Index arrays should have same length"
    assert len(d2d) == len(idx1), "Distance array should match index length"


def test_optimized_no_matches():
    """Test when no matches are found"""
    cata = SkyCoord([0, 1, 2], [0, 1, 2], unit=(u.degree, u.degree))
    grid = SkyCoord([10, 11, 12], [10, 11, 12], unit=(u.degree, u.degree))
    
    r_cut = 1 * u.arcsec
    
    idx1, idx2, _, _ = optimized_search_around_sky(
        cata, grid, r_cut, n_workers=1, verbose=False
    )
    
    assert len(idx1) == 0, "Should find no matches"
    assert len(idx2) == 0, "Should find no matches"


def test_distance_calculations():
    """Test that d2d and d3d are calculated correctly"""
    np.random.seed(789)
    n_cata = 1_000
    n_grid = 500
    
    cata = SkyCoord(
        np.random.uniform(0, 10, n_cata),
        np.random.uniform(0, 10, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(2, 8, n_grid),
        np.random.uniform(2, 8, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 60 * u.arcsec
    
    # Test with compute_3d=False (default)
    idx1, idx2, d2d, d3d = optimized_search_around_sky(
        cata, grid, r_cut, n_workers=2, compute_3d=False, verbose=False
    )
    
    # Verify d2d matches actual separation
    if len(idx1) > 0:
        for i in range(min(10, len(idx1))):  # Check first 10 pairs
            actual_sep = cata[idx1[i]].separation(grid[idx2[i]])
            assert np.abs(d2d[i].arcsec - actual_sep.arcsec) < 0.01, \
                f"d2d mismatch at pair {i}: {d2d[i].arcsec} vs {actual_sep.arcsec}"
        
        # When compute_3d=False, d3d should be 2*sin(d2d/2)
        expected_d3d = 2 * np.sin(0.5 * d2d)
        assert np.allclose(d3d.value, expected_d3d.value, rtol=1e-10), \
            "d3d should be 2*sin(d2d/2) when compute_3d=False"
    
    # Test with compute_3d=True
    idx1b, idx2b, d2d_b, d3d_b = optimized_search_around_sky(
        cata, grid, r_cut, n_workers=2, compute_3d=True, verbose=False
    )
    
    # Results should be the same
    assert len(idx1) == len(idx1b), "compute_3d should not affect matching"
    assert np.array_equal(idx1, idx1b), "Indices should be identical"
    assert np.array_equal(idx2, idx2b), "Indices should be identical"
    assert np.allclose(d2d.arcsec, d2d_b.arcsec), "d2d should be identical"


def test_distance_vs_original():
    """Test that distances match astropy's implementation"""
    np.random.seed(999)
    n_cata = 5_000
    n_grid = 2_000
    
    cata = SkyCoord(
        np.random.uniform(0, 5, n_cata),
        np.random.uniform(0, 5, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(1, 4, n_grid),
        np.random.uniform(1, 4, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 45 * u.arcsec
    
    # Get results from both implementations
    idx1, idx2, d2d, d3d = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, d3d_orig = cata.search_around_sky(grid, r_cut)
    
    # Create dictionaries for easy comparison
    opt_distances = {(i1, i2): d for i1, i2, d in zip(idx1, idx2, d2d)}
    orig_distances = {(i1, i2): d for i1, i2, d in zip(idx1_orig, idx2_orig, d2d_orig)}
    
    # Check that distances match for all pairs
    for pair in opt_distances:
        assert pair in orig_distances, f"Pair {pair} not in original results"
        assert np.abs(opt_distances[pair].arcsec - orig_distances[pair].arcsec) < 0.01, \
            f"Distance mismatch for pair {pair}"


@pytest.mark.parametrize("n_workers", [1, 2, 4, -1])
def test_optimized_different_workers(n_workers):
    """Test that different worker counts produce same results"""
    np.random.seed(456)
    n_cata = 5_000
    n_grid = 2_000
    
    cata = SkyCoord(
        np.random.uniform(0, 5, n_cata),
        np.random.uniform(0, 5, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(1, 4, n_grid),
        np.random.uniform(1, 4, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 45 * u.arcsec
    
    idx1, idx2, d2d, d3d = optimized_search_around_sky(
        cata, grid, r_cut, n_workers=n_workers, verbose=False
    )
    
    # Compare with single worker result
    idx1_ref, idx2_ref, _, _ = optimized_search_around_sky(
        cata, grid, r_cut, n_workers=1, verbose=False
    )
    
    pairs = set(zip(idx1, idx2))
    pairs_ref = set(zip(idx1_ref, idx2_ref))
    
    assert pairs == pairs_ref, f"Results differ with n_workers={n_workers}"


@pytest.mark.parametrize("n1,n2,seed", [
    (1_000, 5_000, 111),   # coords1 smaller
    (5_000, 1_000, 222),   # coords1 larger
    (4_000, 20_000, 333),  # realistic: smaller vs larger
])
def test_kdtree_auto_selection(n1, n2, seed):
    """Test that KDTree is automatically built with smaller array"""
    np.random.seed(seed)
    
    coords1 = SkyCoord(
        np.random.uniform(0, 10, n1),
        np.random.uniform(0, 10, n1),
        unit=(u.degree, u.degree)
    )
    
    coords2 = SkyCoord(
        np.random.uniform(2, 8, n2),
        np.random.uniform(2, 8, n2),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 60 * u.arcsec
    
    # Run optimized version
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        coords1, coords2, r_cut, n_workers=2, verbose=False
    )
    
    # Compare with original astropy
    # coords1.search_around_sky(coords2) internally calls search_around_sky(coords2, coords1)
    # and returns (idx_coords2, idx_coords1)
    idx2_orig, idx1_orig, d2d_orig, _ = coords1.search_around_sky(coords2, r_cut)
    
    # Verify results match
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"Results mismatch: {len(pairs_opt ^ pairs_orig)} pairs differ"
    
    # Verify distances match for all pairs
    if len(idx1) > 0:
        opt_distances = {(i1, i2): d for i1, i2, d in zip(idx1, idx2, d2d)}
        orig_distances = {(i1, i2): d for i1, i2, d in zip(idx1_orig, idx2_orig, d2d_orig)}
        
        for pair in opt_distances:
            assert pair in orig_distances, f"Pair {pair} not in original results"
            assert np.abs(opt_distances[pair].arcsec - orig_distances[pair].arcsec) < 0.01, \
                f"Distance mismatch for pair {pair}"


def test_north_pole_region():
    """Test coordinates near north pole where RA convergence matters"""
    np.random.seed(111)
    n_cata = 2_000
    n_grid = 1_000
    
    # Near north pole (Dec 85-90°)
    cata = SkyCoord(
        np.random.uniform(0, 360, n_cata),
        np.random.uniform(85, 90, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(0, 360, n_grid),
        np.random.uniform(87, 90, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 120 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"North pole: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_south_pole_region():
    """Test coordinates near south pole"""
    np.random.seed(222)
    n_cata = 2_000
    n_grid = 1_000
    
    # Near south pole (Dec -90 to -85°)
    cata = SkyCoord(
        np.random.uniform(0, 360, n_cata),
        np.random.uniform(-90, -85, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(0, 360, n_grid),
        np.random.uniform(-90, -87, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 120 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"South pole: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_ra_wraparound():
    """Test coordinates crossing RA=0°/360° boundary"""
    np.random.seed(333)
    n_cata = 3_000
    n_grid = 1_500
    
    # Coordinates around RA=0° (358° to 2°)
    ra_cata = np.random.uniform(358, 362, n_cata) % 360
    ra_grid = np.random.uniform(358, 362, n_grid) % 360
    
    cata = SkyCoord(
        ra_cata,
        np.random.uniform(20, 40, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        ra_grid,
        np.random.uniform(25, 35, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 90 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"RA wraparound: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_southern_hemisphere():
    """Test coordinates in southern hemisphere"""
    np.random.seed(444)
    n_cata = 5_000
    n_grid = 2_000
    
    cata = SkyCoord(
        np.random.uniform(0, 360, n_cata),
        np.random.uniform(-60, -30, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(0, 360, n_grid),
        np.random.uniform(-50, -35, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 75 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"Southern hemisphere: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_allsky_uniform():
    """Test uniform distribution across entire sky"""
    np.random.seed(555)
    n_cata = 10_000
    n_grid = 5_000
    
    # Uniform sphere distribution using inverse CDF method
    cata = SkyCoord(
        np.random.uniform(0, 360, n_cata),
        np.degrees(np.arcsin(np.random.uniform(-1, 1, n_cata))),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(0, 360, n_grid),
        np.degrees(np.arcsin(np.random.uniform(-1, 1, n_grid))),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 60 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"All-sky: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_clustered_distribution():
    """Test clustered (Gaussian) distribution"""
    np.random.seed(666)
    n_cata = 5_000
    n_grid = 2_000
    
    # Clustered around RA=180°, Dec=30°
    cata = SkyCoord(
        np.random.normal(180, 5, n_cata) % 360,
        np.clip(np.random.normal(30, 2, n_cata), -90, 90),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.normal(180, 3, n_grid) % 360,
        np.clip(np.random.normal(30, 1.5, n_grid), -90, 90),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 90 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"Clustered: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_galactic_coordinates():
    """Test with Galactic coordinate system"""
    np.random.seed(777)
    n_cata = 3_000
    n_grid = 1_500
    
    # Galactic coordinates near plane
    cata = SkyCoord(
        l=np.random.uniform(0, 360, n_cata),
        b=np.random.uniform(-10, 10, n_cata),
        unit=(u.degree, u.degree),
        frame='galactic'
    )
    
    grid = SkyCoord(
        l=np.random.uniform(50, 150, n_grid),
        b=np.random.uniform(-5, 5, n_grid),
        unit=(u.degree, u.degree),
        frame='galactic'
    )
    
    r_cut = 60 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"Galactic coords: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_mixed_coordinate_frames():
    """Test with different coordinate frames (ICRS vs Galactic)"""
    np.random.seed(888)
    n_cata = 2_000
    n_grid = 1_000
    
    # ICRS coordinates
    cata = SkyCoord(
        np.random.uniform(100, 200, n_cata),
        np.random.uniform(20, 40, n_cata),
        unit=(u.degree, u.degree),
        frame='icrs'
    )
    
    # Galactic coordinates (will be transformed)
    grid = SkyCoord(
        l=np.random.uniform(150, 250, n_grid),
        b=np.random.uniform(30, 50, n_grid),
        unit=(u.degree, u.degree),
        frame='galactic'
    )
    
    r_cut = 90 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"Mixed frames: {len(pairs_opt ^ pairs_orig)} pairs differ"


@pytest.mark.parametrize("r_cut", [
    0.1 * u.arcsec,   # Very small
    1 * u.arcsec,     # Small
    30 * u.arcsec,    # Medium-small
    120 * u.arcsec,   # Medium
    300 * u.arcsec,   # Large
    0.5 * u.degree,   # Very large
])
def test_various_search_radii(r_cut):
    """Test with various search radii from tiny to large"""
    np.random.seed(999)
    n_cata = 2_000
    n_grid = 1_000
    
    cata = SkyCoord(
        np.random.uniform(100, 110, n_cata),
        np.random.uniform(30, 40, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(102, 108, n_grid),
        np.random.uniform(32, 38, n_grid),
        unit=(u.degree, u.degree)
    )
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, \
        f"Radius {r_cut}: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_empty_coords():
    """Test with empty coordinate arrays"""
    cata = SkyCoord([], [], unit=(u.degree, u.degree))
    grid = SkyCoord([10, 20], [30, 40], unit=(u.degree, u.degree))
    
    r_cut = 60 * u.arcsec
    
    idx1, idx2, _, _ = optimized_search_around_sky(
        cata, grid, r_cut, n_workers=1, verbose=False
    )
    
    assert len(idx1) == 0, "Empty coords should return no matches"
    assert len(idx2) == 0, "Empty coords should return no matches"


def test_single_point():
    """Test with single coordinate point"""
    cata = SkyCoord([100], [30], unit=(u.degree, u.degree))
    grid = SkyCoord([100.001], [30.001], unit=(u.degree, u.degree))
    
    r_cut = 10 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=1, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, "Single point results differ"


def test_extreme_imbalance():
    """Test with extremely imbalanced array sizes"""
    np.random.seed(1010)
    
    # Very small vs large
    cata = SkyCoord(
        np.random.uniform(100, 110, 50_000),
        np.random.uniform(30, 40, 50_000),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(104, 106, 10),
        np.random.uniform(34, 36, 10),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 90 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"Extreme imbalance: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_equator_region():
    """Test coordinates near celestial equator"""
    np.random.seed(1111)
    n_cata = 5_000
    n_grid = 2_000
    
    # Near equator (Dec -5° to +5°)
    cata = SkyCoord(
        np.random.uniform(0, 360, n_cata),
        np.random.uniform(-5, 5, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(0, 360, n_grid),
        np.random.uniform(-3, 3, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 60 * u.arcsec
    
    idx1, idx2, d2d, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, verbose=False
    )
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_opt = set(zip(idx1, idx2))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_opt == pairs_orig, f"Equator: {len(pairs_opt ^ pairs_orig)} pairs differ"


def test_batch_processing():
    """Test batch processing with data larger than batch_size"""
    np.random.seed(1234)
    # Use data size that requires multiple batches
    n_cata = 150_000  # Will be split into 2 batches with batch_size=100_000
    n_grid = 80_000
    
    cata = SkyCoord(
        np.random.uniform(0, 10, n_cata),
        np.random.uniform(0, 10, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(2, 8, n_grid),
        np.random.uniform(2, 8, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 60 * u.arcsec
    
    # Test with small batch_size to force multiple batches
    idx1_batch, idx2_batch, d2d_batch, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, batch_size=100_000, verbose=False
    )
    
    # Compare with original
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_batch = set(zip(idx1_batch, idx2_batch))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_batch == pairs_orig, f"Batch processing: {len(pairs_batch ^ pairs_orig)} pairs differ"
    assert len(idx1_batch) == len(idx1_orig), "Different number of matches"


def test_batch_processing_large():
    """Test batch processing with larger data requiring 3+ batches"""
    np.random.seed(5678)
    n_cata = 350_000  # Will be split into 4 batches with batch_size=100_000
    n_grid = 150_000
    
    cata = SkyCoord(
        np.random.uniform(0, 10, n_cata),
        np.random.uniform(0, 10, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(2, 8, n_grid),
        np.random.uniform(2, 8, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 45 * u.arcsec
    
    # Test with batch_size=100_000 (4 batches)
    idx1_batch, idx2_batch, d2d_batch, _ = optimized_search_around_sky(
        grid, cata, r_cut, n_workers=2, batch_size=100_000, verbose=False
    )
    
    # Compare with original
    idx1_orig, idx2_orig, d2d_orig, _ = cata.search_around_sky(grid, r_cut)
    
    pairs_batch = set(zip(idx1_batch, idx2_batch))
    pairs_orig = set(zip(idx1_orig, idx2_orig))
    
    assert pairs_batch == pairs_orig, f"Large batch processing: {len(pairs_batch ^ pairs_orig)} pairs differ"
    
    # Verify distances match
    if len(idx1_batch) > 0:
        batch_distances = {(i1, i2): d for i1, i2, d in zip(idx1_batch, idx2_batch, d2d_batch)}
        orig_distances = {(i1, i2): d for i1, i2, d in zip(idx1_orig, idx2_orig, d2d_orig)}
        
        for pair in list(batch_distances.keys())[:100]:  # Check first 100 pairs
            assert pair in orig_distances, f"Pair {pair} not in original results"
            assert np.abs(batch_distances[pair].arcsec - orig_distances[pair].arcsec) < 0.01, \
                f"Distance mismatch for pair {pair}"


if __name__ == "__main__":
    # Run with verbose output for manual testing
    import time
    
    print("Running performance comparison test...")
    print("="*60)
    
    np.random.seed(42)
    n_cata = 100_000
    n_grid = 20_000
    
    cata = SkyCoord(
        np.random.uniform(0, 30, n_cata),
        np.random.uniform(20, 40, n_cata),
        unit=(u.degree, u.degree)
    )
    
    grid = SkyCoord(
        np.random.uniform(10, 20, n_grid),
        np.random.uniform(25, 35, n_grid),
        unit=(u.degree, u.degree)
    )
    
    r_cut = 90 * u.arcsec
    
    print(f"Test data: {n_cata:,} spectra, {n_grid:,} grid points")
    print(f"Search radius: {r_cut}")
    
    # Test optimized version
    print("\nOptimized version:")
    t0 = time.time()
    idx_cata, idx_grid, d2d, d3d = optimized_search_around_sky(cata, grid, r_cut, n_workers=4)
    elapsed = time.time() - t0
    print(f"Time: {elapsed:.2f}s, found {len(idx_cata):,} pairs")
    
    # Test original version
    print("\nOriginal version:")
    t0 = time.time()
    # cata.search_around_sky(grid) internally calls search_around_sky(grid, cata)
    # and returns (idx_grid, idx_cata)
    idx_grid_orig, idx_cata_orig, d2d_orig, d3d_orig = cata.search_around_sky(grid, r_cut)
    elapsed_orig = time.time() - t0
    print(f"Time: {elapsed_orig:.2f}s, found {len(idx_cata_orig):,} pairs")
    
    # Verify results - compare (cata_idx, grid_idx) pairs
    pairs_opt = set(zip(idx_cata, idx_grid))
    pairs_orig = set(zip(idx_cata_orig, idx_grid_orig))
    
    if pairs_opt == pairs_orig:
        print("\n✓ Results verified!")
        print(f"Speedup: {elapsed_orig/elapsed:.2f}x")
    else:
        print(f"\n✗ Results mismatch: {len(pairs_opt ^ pairs_orig)} pairs differ")
