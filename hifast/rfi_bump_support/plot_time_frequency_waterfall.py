#!/usr/bin/env python3
"""Plot raw and processed time-frequency waterfalls for one beam/polarization."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter1d, median_filter

from detect_group_x_bump_sequences import (
    baseline_residual,
    best_sequence,
    detect_bumps,
    odd_window,
    robust_scale,
)


def frequency_slice(
    frequency: np.ndarray,
    minimum: float,
    maximum: float,
) -> slice:
    indices = np.flatnonzero(
        (frequency >= minimum) & (frequency <= maximum)
    )
    if not len(indices):
        raise ValueError(f"No channels in {minimum}–{maximum} MHz")
    return slice(int(indices[0]), int(indices[-1]) + 1)


def fill_rows(values: np.ndarray) -> np.ndarray:
    result = values.astype(np.float64, copy=True)
    channel_indices = np.arange(values.shape[1])
    for row in range(values.shape[0]):
        finite = np.isfinite(result[row])
        if np.count_nonzero(finite) < 2:
            result[row] = np.nan
        elif not np.all(finite):
            result[row, ~finite] = np.interp(
                channel_indices[~finite],
                channel_indices[finite],
                result[row, finite],
            )
    return result


def subtract_local_baseline(
    values: np.ndarray,
    delta_frequency_mhz: float,
    width_mhz: float,
) -> np.ndarray:
    filled = fill_rows(values)
    baseline = median_filter(
        filled,
        size=(1, odd_window(width_mhz, delta_frequency_mhz)),
        mode="nearest",
    )
    residual = values - baseline
    residual[~np.isfinite(values)] = np.nan
    return residual


def bin_time(
    values: np.ndarray,
    mjd: np.ndarray,
    samples_per_bin: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    row_count = (len(values) + samples_per_bin - 1) // samples_per_bin
    binned = np.full((row_count, values.shape[1]), np.nan, dtype=np.float64)
    binned_mjd = np.full(row_count, np.nan, dtype=np.float64)
    bin_start_offsets = np.empty(row_count, dtype=np.int64)
    bin_stop_offsets = np.empty(row_count, dtype=np.int64)
    for row in range(row_count):
        start = row * samples_per_bin
        stop = min((row + 1) * samples_per_bin, len(values))
        bin_start_offsets[row] = start
        bin_stop_offsets[row] = stop
        binned[row] = np.nanmean(values[start:stop], axis=0)
        binned_mjd[row] = np.nanmean(mjd[start:stop])
    return binned, binned_mjd, bin_start_offsets, bin_stop_offsets


def symmetric_limit(values: np.ndarray, percentile: float = 99.0) -> float:
    finite = np.abs(values[np.isfinite(values)])
    return float(np.percentile(finite, percentile))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-file", type=Path, required=True)
    parser.add_argument("--final-file", type=Path, required=True)
    parser.add_argument("--cube", type=Path, required=True)
    parser.add_argument("--polar", choices=["XX", "YY"], default="YY")
    parser.add_argument("--x-range", type=float, nargs=2, required=True)
    parser.add_argument(
        "--reference-x-range", type=float, nargs=2, required=True
    )
    parser.add_argument(
        "--freq-range", type=float, nargs=2, default=[1421.0, 1424.8]
    )
    parser.add_argument("--time-bin-seconds", type=float, default=15.0)
    parser.add_argument("--baseline-width-mhz", type=float, default=0.55)
    parser.add_argument("--smooth-fwhm-mhz", type=float, default=0.07)
    parser.add_argument("--output-npz", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()

    celestial_wcs = WCS(fits.getheader(args.cube, 0)).celestial
    with h5py.File(args.final_file, "r") as final_handle:
        ra = final_handle["S/ra"][()]
        dec = final_handle["S/dec"][()]
        final_mjd = final_handle["S/mjd"][()]
    x_pixel, _ = celestial_wcs.world_to_pixel_values(ra, dec)

    selected = (
        (x_pixel >= args.x_range[0]) & (x_pixel <= args.x_range[1])
    )
    reference = (
        (x_pixel >= args.reference_x_range[0])
        & (x_pixel <= args.reference_x_range[1])
    )
    selected_indices = np.flatnonzero(selected)
    reference_indices = np.flatnonzero(reference)
    if not len(selected_indices) or not len(reference_indices):
        raise RuntimeError("The requested x range has no time samples")

    with h5py.File(args.stage_file, "r") as stage_handle:
        spectra = stage_handle["S"]
        frequency_all = spectra["freq"][()]
        channel_slice = frequency_slice(
            frequency_all, args.freq_range[0], args.freq_range[1]
        )
        frequency = frequency_all[channel_slice]
        polarization = 0 if args.polar == "XX" else 1
        values = spectra["Ta"][
            polarization, selected_indices, channel_slice
        ].astype(np.float64)
        reference_values = spectra["Ta"][
            polarization, reference_indices, channel_slice
        ].astype(np.float64)
        if "is_rfi" in spectra:
            rfi = spectra["is_rfi"][selected_indices, channel_slice]
            reference_rfi = spectra["is_rfi"][
                reference_indices, channel_slice
            ]
            values[rfi] = np.nan
            reference_values[reference_rfi] = np.nan

    mjd = final_mjd[selected_indices]
    reference_mjd = final_mjd[reference_indices]
    delta_frequency = float(np.median(np.diff(frequency)))
    sample_seconds = float(np.median(np.diff(mjd)) * 86400.0)
    samples_per_bin = max(
        1, int(round(args.time_bin_seconds / sample_seconds))
    )

    raw_residual = subtract_local_baseline(
        values, delta_frequency, args.baseline_width_mhz
    )
    (
        binned_values,
        binned_mjd,
        bin_start_offsets,
        bin_stop_offsets,
    ) = bin_time(
        values, mjd, samples_per_bin
    )
    binned_residual = subtract_local_baseline(
        binned_values, delta_frequency, args.baseline_width_mhz
    )
    sigma_channel = (
        args.smooth_fwhm_mhz / 2.354820045 / abs(delta_frequency)
    )
    smoothed = gaussian_filter1d(
        fill_rows(binned_residual),
        sigma=sigma_channel,
        axis=1,
        mode="nearest",
    )
    smoothed_scale = robust_scale(smoothed)
    smoothed_standardized = smoothed / smoothed_scale

    full_mean = np.nanmean(values, axis=0)
    reference_mean = np.nanmean(reference_values, axis=0)
    full_mean_residual = baseline_residual(
        full_mean, delta_frequency, "local_median"
    )
    reference_mean_residual = baseline_residual(
        reference_mean, delta_frequency, "local_median"
    )
    chain = best_sequence(
        detect_bumps(frequency, reference_mean_residual)
    )
    blind_candidate_centers = np.asarray(
        [bump.frequency_mhz for bump in chain], dtype=np.float64
    )
    blind_candidate_widths = np.asarray(
        [bump.fwhm_mhz for bump in chain], dtype=np.float64
    )
    candidate_centers = blind_candidate_centers.copy()
    if len(blind_candidate_centers) >= 2:
        comb_spacing = float(np.median(np.diff(blind_candidate_centers)))
        reference_smoothed = gaussian_filter1d(
            reference_mean_residual,
            sigma=sigma_channel,
            mode="nearest",
        )
        expected_centers = []
        anchor = float(blind_candidate_centers[0])
        for offset in range(-10, 11):
            expected = anchor + offset * comb_spacing
            if not (
                frequency[0] + 0.30
                <= expected
                <= frequency[-1] - 0.30
            ):
                continue
            search = np.flatnonzero(np.abs(frequency - expected) <= 0.08)
            channel = int(search[np.nanargmax(reference_smoothed[search])])
            expected_centers.append(float(frequency[channel]))
        candidate_centers = np.asarray(
            sorted(set(expected_centers)), dtype=np.float64
        )

    elapsed_minutes = (mjd - mjd[0]) * 1440.0
    binned_minutes = (binned_mjd - mjd[0]) * 1440.0
    reference_minutes = (
        np.asarray([reference_mjd[0], reference_mjd[-1]]) - mjd[0]
    ) * 1440.0

    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        frequency_mhz=frequency,
        selected_original_indices=selected_indices,
        selected_mjd=mjd,
        elapsed_minutes=elapsed_minutes,
        binned_mjd=binned_mjd,
        binned_elapsed_minutes=binned_minutes,
        bin_start_offsets=bin_start_offsets,
        bin_stop_offsets=bin_stop_offsets,
        raw_residual_k=raw_residual.astype(np.float32),
        binned_residual_k=binned_residual.astype(np.float32),
        smoothed_standardized=smoothed_standardized.astype(np.float32),
        full_mean_residual_k=full_mean_residual.astype(np.float32),
        reference_mean_residual_k=reference_mean_residual.astype(np.float32),
        candidate_centers_mhz=candidate_centers,
        blind_candidate_centers_mhz=blind_candidate_centers,
        blind_candidate_widths_mhz=blind_candidate_widths,
        x_range=np.asarray(args.x_range),
        reference_x_range=np.asarray(args.reference_x_range),
        sample_seconds=sample_seconds,
        samples_per_bin=samples_per_bin,
        polarization=args.polar,
    )

    figure, axes = plt.subplots(
        4,
        1,
        figsize=(15, 13),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.15, 1.15, 1.15]},
    )
    axes[0].plot(
        frequency,
        full_mean_residual * 1000,
        color="0.45",
        lw=0.8,
        label=f"full x={args.x_range[0]:g}–{args.x_range[1]:g}",
    )
    axes[0].plot(
        frequency,
        reference_mean_residual * 1000,
        color="tab:red",
        lw=1.1,
        label=(
            f"reference x={args.reference_x_range[0]:g}–"
            f"{args.reference_x_range[1]:g}"
        ),
    )
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_ylabel("mean residual\n(mK)")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.2)

    raw_limit = symmetric_limit(raw_residual * 1000, 99.0)
    raw_image = axes[1].imshow(
        raw_residual * 1000,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=[
            frequency[0],
            frequency[-1],
            elapsed_minutes[0],
            elapsed_minutes[-1],
        ],
        cmap="RdBu_r",
        vmin=-raw_limit,
        vmax=raw_limit,
    )
    axes[1].set_ylabel("time from start\n(min)")
    axes[1].set_title(
        f"Every {sample_seconds:.3f} s spectrum after "
        f"{args.baseline_width_mhz:.2f} MHz local-baseline subtraction"
    )
    figure.colorbar(raw_image, ax=axes[1], label="residual (mK)", pad=0.01)

    binned_limit = symmetric_limit(binned_residual * 1000, 99.0)
    binned_image = axes[2].imshow(
        binned_residual * 1000,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=[
            frequency[0],
            frequency[-1],
            binned_minutes[0],
            binned_minutes[-1],
        ],
        cmap="RdBu_r",
        vmin=-binned_limit,
        vmax=binned_limit,
    )
    axes[2].set_ylabel("time from start\n(min)")
    axes[2].set_title(
        f"Time-binned residual; {samples_per_bin} samples "
        f"≈{samples_per_bin * sample_seconds:.1f} s per row"
    )
    figure.colorbar(
        binned_image, ax=axes[2], label="residual (mK)", pad=0.01
    )

    standardized_limit = max(
        2.0, symmetric_limit(smoothed_standardized, 99.0)
    )
    smooth_image = axes[3].imshow(
        smoothed_standardized,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=[
            frequency[0],
            frequency[-1],
            binned_minutes[0],
            binned_minutes[-1],
        ],
        cmap="RdBu_r",
        vmin=-standardized_limit,
        vmax=standardized_limit,
    )
    axes[3].set_ylabel("time from start\n(min)")
    axes[3].set_xlabel("topocentric frequency (MHz)")
    axes[3].set_title(
        f"Same time bins plus {args.smooth_fwhm_mhz:.3f} MHz "
        "frequency smoothing"
    )
    figure.colorbar(
        smooth_image,
        ax=axes[3],
        label="smoothed response / robust RMS",
        pad=0.01,
    )

    for axis in axes:
        for center_index, center in enumerate(candidate_centers):
            axis.axvline(
                center,
                color="tab:red",
                lw=0.8,
                linestyle="--",
                alpha=0.8,
                label=(
                    "comb guide from measured spacing"
                    if axis is axes[0] and center_index == 0
                    else None
                ),
            )
    for axis in axes[1:]:
        axis.axhspan(
            reference_minutes[0],
            reference_minutes[1],
            color="gold",
            alpha=0.10,
        )
    axes[0].legend(fontsize=8, ncol=3)

    figure.suptitle(
        f"Scan 43, M05 {args.polar}, after standing-wave removal; "
        f"cube x={args.x_range[0]:g}–{args.x_range[1]:g}",
        y=0.995,
    )
    figure.tight_layout(rect=[0, 0, 1, 0.985])
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=170)
    plt.close(figure)

    centers_text = ";".join(f"{value:.9f}" for value in candidate_centers)
    print(f"selected_samples={len(values)}")
    print(f"duration_minutes={elapsed_minutes[-1]:.6f}")
    print(f"samples_per_bin={samples_per_bin}")
    print(f"candidate_centers_mhz={centers_text}")
    print(f"saved_npz={args.output_npz}")
    print(f"saved_figure={args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
