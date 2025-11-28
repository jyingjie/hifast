"""
Optimized fill_grid implementation

Uses query_ball_point instead of query_ball_tree for significant performance 
and memory improvements.

Key advantages:
1. Only builds one KDTree (grid) instead of two
2. Leverages scipy's built-in multi-threading (workers parameter)
3. 63% memory reduction (20 GB vs 54.5 GB)
4. 6.4x performance improvement with 10 workers (52s vs 330s)
5. Uses threadpool_limits to control BLAS/OpenMP threads consistently

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
from threadpoolctl import threadpool_limits
import gc


def optimized_search_around_sky(coords1, coords2, seplimit, n_workers=1, 
                                compute_3d=False, compute_2d=True, verbose=True, 
                                batch_size=1_000_000, sep_batch_size=10_000_000):
    """
    Optimized search_around_sky implementation with automatic KDTree size optimization.
    
    Uses query_ball_point with batch processing to minimize memory usage.
    Automatically builds KDTree with the smaller coordinate set.
    Uses threadpool_limits to control all BLAS/OpenMP threads consistently with n_workers.
    
    Parameters
    ----------
    coords1 : SkyCoord
        First set of coordinates
    coords2 : SkyCoord
        Second set of coordinates
    seplimit : Quantity
        Angular separation limit, must be scalar
    n_workers : int, optional
        Number of parallel threads. Default is 1 (single-threaded).
        Set to -1 to use all CPU cores.
        This controls both scipy workers and BLAS/OpenMP threads.
    compute_3d : bool, optional
        Whether to compute 3D distances. Default is False to save time
    compute_2d : bool, optional
        Whether to compute 2D angular separations. Default is True.
        Set to False to skip separation calculation and save memory/time
        if you only need the match indices.
    verbose : bool, optional
        Whether to print progress messages. Default is True
    batch_size : int, optional
        Number of query points to process in each batch. Default is 1,000,000.
        Larger values are faster but use more memory.
        Recommended: 500,000 (low memory) to 2,000,000 (high memory).
    sep_batch_size : int, optional
        Number of matches to process in each batch when computing separations.
        Default is 10,000,000. Only used if compute_2d=True.
        Larger values are faster but use more memory during separation calculation.
    
    Returns
    -------
    idx1 : ndarray
        Indices into coords1, sorted in ascending order
    idx2 : ndarray
        Indices into coords2, corresponding to idx1
    sep2d : Angle
        2D angular separation for each pair
    dist3d : Quantity
        3D distance (only if compute_3d=True, otherwise returns 2*sin(sep2d/2))
    
    Notes
    -----
    - Results are sorted by idx1 (indices into coords1) in ascending order
    - Automatically builds KDTree with the smaller coordinate set to save memory
    - Uses batch processing to reduce peak memory usage
    - Only supports scalar seplimit
    - Uses scipy's multi-threading, more efficient than multiprocessing
    - Uses threadpool_limits to ensure consistent thread usage across all operations
    - Return format matches astropy's search_around_sky
    
    Examples
    --------
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> coords_a = SkyCoord(ra_array_a, dec_array_a, unit=(u.degree, u.degree))
    >>> coords_b = SkyCoord(ra_array_b, dec_array_b, unit=(u.degree, u.degree))
    >>> # Use 10 threads with default batch size
    >>> idx_a, idx_b, d2d, d3d = optimized_search_around_sky(
    ...     coords_a, coords_b, 90*u.arcsec, n_workers=10
    ... )
    >>> # Use smaller batch size for low memory systems
    >>> idx_a, idx_b, d2d, d3d = optimized_search_around_sky(
    ...     coords_a, coords_b, 90*u.arcsec, n_workers=10, batch_size=500_000
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
    
    # Determine number of parallel threads
    if n_workers < 0:
        n_workers = cpu_count()
    elif n_workers == 0:
        n_workers = 1
    
    # Use threadpool_limits to control all thread-based parallelism
    # This ensures numpy, scipy, and other libraries use consistent thread count
    # user_api=None (default) applies to both BLAS and OpenMP libraries
    with threadpool_limits(limits=n_workers):
        return _optimized_search_around_sky_impl(
            coords1, coords2, seplimit, n_workers, compute_3d, compute_2d, 
            verbose, batch_size, sep_batch_size
        )


def _optimized_search_around_sky_impl(coords1, coords2, seplimit, n_workers, 
                                      compute_3d, compute_2d, verbose, batch_size,
                                      sep_batch_size):
    """
    Internal implementation of optimized_search_around_sky with batch processing.
    This function runs within threadpool_limits context.
    """
    # Transform coordinate frames first (always transform coords1 to coords2's frame)
    coords1 = coords1.transform_to(coords2)
    
    # Determine which array is smaller to build KDTree (saves memory)
    # Swap if needed so that tree_coords is always the smaller one
    if len(coords1) > len(coords2):
        # coords2 is smaller, swap them
        tree_coords = coords2
        query_coords = coords1
        swapped = True
    else:
        # coords1 is smaller or equal, use as-is
        tree_coords = coords1
        query_coords = coords2
        swapped = False
    
    # Convert to Cartesian coordinates sequentially to reduce peak memory
    # Process tree_coords first (smaller array for KDTree)
    urepr_tree = tree_coords.data.represent_as(UnitSphericalRepresentation)
    tree_coords_frame = tree_coords.realize_frame(urepr_tree)
    cart_tree = tree_coords_frame.cartesian
    cartxyz_tree = np.array([cart_tree.x.value, cart_tree.y.value, cart_tree.z.value]).T
    
    # Clean up tree intermediate objects immediately
    del urepr_tree, tree_coords_frame, cart_tree
    
    # Now process query_coords (larger array)
    urepr_query = query_coords.data.represent_as(UnitSphericalRepresentation)
    query_coords_frame = query_coords.realize_frame(urepr_query)
    cart_query = query_coords_frame.cartesian
    cartxyz_query = np.array([cart_query.x.value, cart_query.y.value, cart_query.z.value]).T
    
    # Clean up query intermediate objects immediately
    del urepr_query, query_coords_frame, cart_query
    
    # Calculate search radius in Cartesian space
    r = (2 * np.sin(Angle(0.5 * seplimit))).value
    
    # Calculate number of batches
    n_query_points = len(query_coords)
    n_batches = (n_query_points + batch_size - 1) // batch_size
    
    if verbose:
        print(f"Optimized fill_grid: using {n_workers} threads (BLAS/OpenMP limited)")
        print(f"  coords1: {len(coords1):,}, coords2: {len(coords2):,}, "
              f"Radius: {seplimit}")
        if swapped:
            print(f"  Building KDTree with coords2 ({len(coords2):,} points, smaller)")
        else:
            print(f"  Building KDTree with coords1 ({len(coords1):,} points, smaller)")
        print(f"  Batch processing: {n_batches} batches of {batch_size:,} points")
    
    # Build KDTree with smaller array (tree_coords)
    kdt = cKDTree(cartxyz_tree)
    
    # Clean up tree Cartesian array after building KDTree
    del cartxyz_tree
    
    # Determine optimal integer dtype based on array sizes
    max_idx = max(len(coords1), len(coords2))
    if max_idx < 2_147_483_647:
        idx_dtype = np.int32  # Save 50% memory for index arrays
    else:
        idx_dtype = np.int64  # Use int64 for very large arrays
    
    # Process in batches to reduce peak memory
    all_idxs_query = []
    all_idxs_tree = []
    
    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, n_query_points)
        batch_len = end_idx - start_idx
        
        if verbose and n_batches > 1:
            print(f"  Processing batch {batch_idx + 1}/{n_batches} "
                  f"({start_idx:,} to {end_idx:,})...")
        
        # Query this batch
        cartxyz_batch = cartxyz_query[start_idx:end_idx]
        matches_batch = kdt.query_ball_point(
            cartxyz_batch,
            r,
            workers=n_workers,
            return_sorted=False  # Don't sort, we'll sort later by coords1 indices
        )
        
        # Immediately convert to index arrays and release matches_batch
        batch_match_counts = np.array([len(m) for m in matches_batch], dtype=np.int32)
        total_matches_batch = np.sum(batch_match_counts, dtype=np.int64)
        
        if total_matches_batch > 0:
            # Build query indices (with offset for batch start)
            idxs_query_batch = np.repeat(
                np.arange(start_idx, end_idx, dtype=idx_dtype),
                batch_match_counts
            )
            
            # Build tree indices
            idxs_tree_batch = np.concatenate(
                [m for m in matches_batch if len(m) > 0]
            ).astype(idx_dtype, copy=False)
            
            all_idxs_query.append(idxs_query_batch)
            all_idxs_tree.append(idxs_tree_batch)
        
        # Clean up batch data immediately
        del matches_batch, batch_match_counts, cartxyz_batch
        if total_matches_batch > 0:
            del idxs_query_batch, idxs_tree_batch
    gc.collect()
    
    # Clean up remaining data
    del cartxyz_query, kdt
    gc.collect()
    
    # Combine all batches
    if len(all_idxs_query) == 0:
        # No matches found
        return (np.array([], dtype=idx_dtype), 
                np.array([], dtype=idx_dtype),
                Angle([], unit=u.degree),
                u.Quantity([], unit=u.dimensionless_unscaled))
    
    idxs_query = np.concatenate(all_idxs_query)
    # Immediately release all_idxs_query to free memory before next concatenate
    del all_idxs_query
    gc.collect()
    
    idxs_tree = np.concatenate(all_idxs_tree)
    # Release all_idxs_tree
    del all_idxs_tree
    gc.collect()
    
    # Report total matches
    total_matches = len(idxs_query)
    if verbose:
        print(f"  Found {total_matches:,} matching pairs")
    
    # Map back to original coords1/coords2 order
    if swapped:
        # We swapped: tree=coords2, query=coords1
        # So idxs_tree are coords2 indices, idxs_query are coords1 indices
        idxs1 = idxs_query
        idxs2 = idxs_tree
    else:
        # No swap: tree=coords1, query=coords2
        # So idxs_tree are coords1 indices, idxs_query are coords2 indices
        idxs1 = idxs_tree
        idxs2 = idxs_query
    
    # Clean up temporary index arrays after mapping
    del idxs_tree, idxs_query
    
    # Sort by coords1 indices (consistent with astropy behavior)
    # Note: indexing preserves dtype, so no need for explicit astype
    sort_order = np.argsort(idxs1)
    idxs1 = idxs1[sort_order]
    idxs2 = idxs2[sort_order]
    
    # Clean up sort order array
    del sort_order
    gc.collect()
    
    # Calculate angular separation and 3D distance
    if not compute_2d:
        # Skip separation calculation, return empty Angle arrays
        if verbose:
            print("  Skipping separation calculation (compute_2d=False)")
            print("  Done!")
        return (idxs1, idxs2,
                Angle([], unit=u.degree),
                u.Quantity([], unit=u.dimensionless_unscaled))
    
    # Compute separations in batches to reduce memory usage
    n_sep_batches = (total_matches + sep_batch_size - 1) // sep_batch_size
    
    if verbose and n_sep_batches > 1:
        print(f"  Computing separations in {n_sep_batches} batches...")
    
    d2d_list = []
    d3d_list = [] if compute_3d else None
    
    for i in range(0, total_matches, sep_batch_size):
        end = min(i + sep_batch_size, total_matches)
        
        if verbose and n_sep_batches > 1:
            batch_num = i // sep_batch_size + 1
            print(f"    Separation batch {batch_num}/{n_sep_batches}...")
        
        # Compute 2D separation for this batch
        d2d_batch = coords1[idxs1[i:end]].separation(coords2[idxs2[i:end]])
        d2d_list.append(d2d_batch.value)
        
        # Compute 3D separation if requested
        if compute_3d:
            try:
                d3d_batch = coords1[idxs1[i:end]].separation_3d(coords2[idxs2[i:end]])
                d3d_list.append(d3d_batch.value)
            except ValueError:
                # No distance info, use unit sphere distance
                d3d_batch = 2 * np.sin(0.5 * d2d_batch)
                d3d_list.append(d3d_batch.value)
        
        # Clean up batch data
        del d2d_batch
        if compute_3d:
            del d3d_batch
    gc.collect()
    
    # Combine all batches
    d2ds = Angle(np.concatenate(d2d_list), unit=u.degree)
    del d2d_list
    gc.collect()
    
    if compute_3d:
        d3ds = u.Quantity(np.concatenate(d3d_list), unit=u.dimensionless_unscaled)
        del d3d_list
    else:
        # Skip 3D calculation to save time
        d3ds = 2 * np.sin(0.5 * d2ds)
    
    gc.collect()
    
    if verbose:
        print("  Done!")
    
    return idxs1, idxs2, d2ds, d3ds


# Backward compatibility alias
parallel_search_around_sky = optimized_search_around_sky
