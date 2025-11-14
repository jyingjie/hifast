"""
Optimized fill_grid implementation

Uses query_ball_point instead of query_ball_tree for significant performance 
and memory improvements.

Key advantages:
1. Only builds one KDTree (grid) instead of two
2. Leverages scipy's built-in multi-threading (workers parameter)
3. 63% memory reduction (20 GB vs 54.5 GB)
4. 6.4x performance improvement with 10 workers (52s vs 330s)

Benchmark (19.7M spectra vs 4M grid points, n_workers=10):
- Original search_around_sky: 330s, peak memory 54.5 GB
- Optimized query_ball_point: 52s, peak memory 20 GB

IMPORTANT NOTE about astropy's search_around_sky:
--------------------------------------------------
The SkyCoord.search_around_sky method has a confusing parameter order:
    
    coords1.search_around_sky(coords2, seplimit)
    
internally calls:
    
    search_around_sky(coords2, coords1, seplimit)  # Note: reversed!
    
and returns (idx_coords2, idx_coords1).

So when replacing cata.search_around_sky(grid) with optimized_search_around_sky,
you must call optimized_search_around_sky(grid, cata) to get the same result order.
"""

import numpy as np
from astropy.coordinates import SkyCoord, Angle
from astropy.coordinates.representation import UnitSphericalRepresentation
from astropy import units as u
from scipy.spatial import cKDTree
from multiprocessing import cpu_count


def optimized_search_around_sky(coords1, coords2, seplimit, n_workers=1, 
                                compute_3d=False, verbose=True):
    """
    Optimized search_around_sky implementation
    
    Uses query_ball_point instead of the original query_ball_tree method.
    
    Parameters
    ----------
    coords1 : SkyCoord
        First set of coordinates (spectra/catalog coordinates)
    coords2 : SkyCoord
        Second set of coordinates (grid coordinates)
    seplimit : Quantity
        Angular separation limit, must be scalar
    n_workers : int, optional
        Number of parallel threads. Default is 1 (single-threaded).
        Set to -1 to use all CPU cores
    compute_3d : bool, optional
        Whether to compute 3D distances. Default is False to save time
    verbose : bool, optional
        Whether to print progress messages. Default is True
    
    Returns
    -------
    idx1 : ndarray
        Indices into coords1 (first coordinate set)
    idx2 : ndarray
        Indices into coords2 (second coordinate set)
    sep2d : Angle
        2D angular separation
    dist3d : Quantity
        3D distance (only if compute_3d=True, otherwise returns 2*sin(sep2d/2))
    
    Notes
    -----
    - Return value order matches astropy's search_around_sky: (idx1, idx2, sep2d, dist3d)
    - Only supports scalar seplimit
    - Uses scipy's multi-threading, more efficient than multiprocessing
    
    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> cata = SkyCoord(ra_array, dec_array, unit=(u.degree, u.degree))
    >>> grid = SkyCoord(grid_ra, grid_dec, unit=(u.degree, u.degree))
    >>> # Use 10 threads
    >>> idx_g, idx_c, d2d, d3d = optimized_search_around_sky(
    ...     cata, grid, 90*u.arcsec, n_workers=10
    ... )
    >>> # Use all CPU cores
    >>> idx_g, idx_c, d2d, d3d = optimized_search_around_sky(
    ...     cata, grid, 90*u.arcsec, n_workers=-1
    ... )
    """
    
    # Parameter validation
    if coords1.ndim != 1 or coords2.ndim != 1:
        raise ValueError("Only supports 1-dimensional coordinate arrays")
    
    if not seplimit.isscalar:
        raise NotImplementedError(
            "Non-scalar seplimit not yet supported. "
            "Use the standard search_around_sky for this case."
        )
    
    # Transform coordinate frames
    coords1 = coords1.transform_to(coords2)
    
    # Convert to unit spherical representation (remove distance info)
    urepr1 = coords1.data.represent_as(UnitSphericalRepresentation)
    urepr2 = coords2.data.represent_as(UnitSphericalRepresentation)
    
    coords1_frame = coords1.realize_frame(urepr1)
    coords2_frame = coords2.realize_frame(urepr2)
    
    # Get Cartesian coordinates
    cartxyz1 = coords1_frame.cartesian.xyz.value.T  # (N, 3)
    cartxyz2 = coords2_frame.cartesian.xyz.value.T  # (M, 3)
    
    # Calculate search radius in Cartesian space
    r = (2 * np.sin(Angle(0.5 * seplimit))).value
    
    # Determine number of parallel threads
    if n_workers < 0:
        n_workers = cpu_count()
    elif n_workers == 0:
        n_workers = 1
    
    if verbose:
        print(f"Optimized fill_grid: using {n_workers} threads")
        print(f"  Spectra: {len(coords1):,}, Grid: {len(coords2):,}, "
              f"Radius: {seplimit}")
    
    # Build KDTree for coords2 (grid)
    # Build the smaller KDTree to save memory
    kdt_coords2 = cKDTree(cartxyz2)
    
    # Query all coords1 points using query_ball_point
    # Leverages scipy's built-in multi-threading
    matches = kdt_coords2.query_ball_point(
        cartxyz1,
        r,
        workers=n_workers  # Multi-threaded parallelism
    )
    
    # Calculate total number of matches
    total_matches = sum(len(m) for m in matches)
    
    if verbose:
        print(f"  Found {total_matches:,} matching pairs")
    
    if total_matches == 0:
        # No matches found
        return (np.array([], dtype=int), 
                np.array([], dtype=int),
                Angle([], unit=u.degree),
                u.Quantity([], unit=u.dimensionless_unscaled))
    
    # Pre-allocate arrays
    idxs1 = np.empty(total_matches, dtype=np.int64)  # coords1 indices
    idxs2 = np.empty(total_matches, dtype=np.int64)  # coords2 indices
    
    # Fill arrays
    pos = 0
    for i, coords2_matches in enumerate(matches):
        n = len(coords2_matches)
        if n > 0:
            idxs1[pos:pos+n] = i                # i is index into coords1
            idxs2[pos:pos+n] = coords2_matches  # matches are indices into coords2
            pos += n
    
    # Sort by coords1 indices (consistent with astropy behavior)
    sort_order = np.argsort(idxs1)
    idxs1 = idxs1[sort_order]
    idxs2 = idxs2[sort_order]
    
    # Calculate angular separation and 3D distance
    d2ds = coords1[idxs1].separation(coords2[idxs2])
    
    if compute_3d:
        try:
            d3ds = coords1[idxs1].separation_3d(coords2[idxs2])
        except ValueError:
            # No distance info, use unit sphere distance
            d3ds = 2 * np.sin(0.5 * d2ds)
    else:
        # Skip 3D calculation to save time
        d3ds = 2 * np.sin(0.5 * d2ds)
    
    if verbose:
        print("  Done!")
    
    return idxs1, idxs2, d2ds, d3ds


# Backward compatibility alias
parallel_search_around_sky = optimized_search_around_sky
