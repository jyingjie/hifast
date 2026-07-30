#!/usr/bin/env python3
"""Build a mapped floating-point Waterfall from a pre-Doppler HDF5 file."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

from detect_group_x_bump_sequences import robust_scale
from plot_time_frequency_waterfall import (
    bin_time,
    fill_rows,
    frequency_slice,
    subtract_local_baseline,
    symmetric_limit,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--field", default="flux")
    parser.add_argument("--polar", choices=["XX", "YY"], required=True)
    parser.add_argument("--time-index-start", type=int, required=True)
    parser.add_argument("--time-index-stop", type=int, required=True)
    parser.add_argument(
        "--freq-range", type=float, nargs=2, default=[1421.0, 1424.8]
    )
    parser.add_argument("--time-bin-seconds", type=float, default=15.0)
    parser.add_argument("--baseline-width-mhz", type=float, default=0.55)
    parser.add_argument("--smooth-fwhm-mhz", type=float, default=0.07)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-npz", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()
    data_unit = "Jy/beam" if args.field == "flux" else "K"
    display_unit = "mJy/beam" if args.field == "flux" else "mK"

    if args.time_index_stop <= args.time_index_start:
        raise ValueError("time-index-stop must be larger than the start")
    original_indices = np.arange(
        args.time_index_start, args.time_index_stop, dtype=np.int64
    )

    with h5py.File(args.input, "r") as handle:
        spectra = handle["S"]
        frequency_all = spectra["freq"][()]
        channel_slice = frequency_slice(
            frequency_all, args.freq_range[0], args.freq_range[1]
        )
        frequency_mhz = frequency_all[channel_slice].astype(np.float64)
        if args.field not in spectra:
            raise KeyError(f"S/{args.field} is not present in {args.input}")
        polarization = 0 if args.polar == "XX" else 1
        values = spectra[args.field][
            polarization,
            args.time_index_start : args.time_index_stop,
            channel_slice,
        ].astype(np.float64)
        mjd = spectra["mjd"][
            args.time_index_start : args.time_index_stop
        ].astype(np.float64)
        existing_rfi = (
            spectra["is_rfi"][
                args.time_index_start : args.time_index_stop,
                channel_slice,
            ].astype(bool)
            if "is_rfi" in spectra
            else np.zeros_like(values, dtype=bool)
        )
        existing_excluded = (
            spectra["is_excluded"][
                args.time_index_start : args.time_index_stop,
                channel_slice,
            ].astype(bool)
            if "is_excluded" in spectra
            else np.zeros_like(values, dtype=bool)
        )

    if len(values) != len(original_indices):
        raise RuntimeError("The HDF5 selection is shorter than requested")
    values[existing_rfi] = np.nan
    delta_frequency_mhz = float(np.median(np.diff(frequency_mhz)))
    sample_seconds = float(np.median(np.diff(mjd)) * 86400.0)
    samples_per_bin = max(
        1, int(round(args.time_bin_seconds / sample_seconds))
    )

    raw_residual = subtract_local_baseline(
        values, delta_frequency_mhz, args.baseline_width_mhz
    )
    (
        binned_values,
        binned_mjd,
        bin_start_offsets,
        bin_stop_offsets,
    ) = bin_time(values, mjd, samples_per_bin)
    binned_residual = subtract_local_baseline(
        binned_values, delta_frequency_mhz, args.baseline_width_mhz
    )
    sigma_channel = (
        args.smooth_fwhm_mhz
        / 2.354820045
        / abs(delta_frequency_mhz)
    )
    smoothed = gaussian_filter1d(
        fill_rows(binned_residual),
        sigma=sigma_channel,
        axis=1,
        mode="nearest",
    )
    smoothed_standardized = smoothed / robust_scale(smoothed)
    elapsed_minutes = (mjd - mjd[0]) * 1440
    binned_elapsed_minutes = (binned_mjd - mjd[0]) * 1440

    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        sample_id=args.sample_id,
        source_hdf5=str(args.input),
        processing_stage="pre_fc_topocentric",
        data_field=args.field,
        polarization=args.polar,
        frequency_mhz=frequency_mhz,
        selected_original_indices=original_indices,
        selected_mjd=mjd,
        elapsed_minutes=elapsed_minutes,
        binned_mjd=binned_mjd,
        binned_elapsed_minutes=binned_elapsed_minutes,
        bin_start_offsets=bin_start_offsets,
        bin_stop_offsets=bin_stop_offsets,
        data_unit=data_unit,
        raw_residual=raw_residual.astype(np.float32),
        binned_residual=binned_residual.astype(np.float32),
        smoothed_standardized=smoothed_standardized.astype(np.float32),
        existing_rfi=existing_rfi,
        existing_excluded=existing_excluded,
        existing_rfi_fraction=np.mean(existing_rfi),
        existing_excluded_fraction=np.mean(existing_excluded),
        sample_seconds=sample_seconds,
        samples_per_bin=samples_per_bin,
        baseline_width_mhz=args.baseline_width_mhz,
        smooth_fwhm_mhz=args.smooth_fwhm_mhz,
    )

    figure, axes = plt.subplots(
        3,
        1,
        figsize=(15, 11),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0]},
    )
    raw_limit = symmetric_limit(raw_residual * 1000, 99.0)
    raw_image = axes[0].imshow(
        raw_residual * 1000,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=[
            frequency_mhz[0],
            frequency_mhz[-1],
            elapsed_minutes[0],
            elapsed_minutes[-1],
        ],
        cmap="RdBu_r",
        vmin=-raw_limit,
        vmax=raw_limit,
    )
    axes[0].set_ylabel("time (min)")
    axes[0].set_title(f"Every {sample_seconds:.3f} s spectrum")
    figure.colorbar(
        raw_image, ax=axes[0], pad=0.01, label=f"residual ({display_unit})"
    )

    binned_limit = symmetric_limit(binned_residual * 1000, 99.0)
    binned_image = axes[1].imshow(
        binned_residual * 1000,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=[
            frequency_mhz[0],
            frequency_mhz[-1],
            binned_elapsed_minutes[0],
            binned_elapsed_minutes[-1],
        ],
        cmap="RdBu_r",
        vmin=-binned_limit,
        vmax=binned_limit,
    )
    axes[1].set_ylabel("time (min)")
    axes[1].set_title(
        f"{samples_per_bin} samples = "
        f"{samples_per_bin * sample_seconds:.1f} s per row"
    )
    figure.colorbar(
        binned_image,
        ax=axes[1],
        pad=0.01,
        label=f"residual ({display_unit})",
    )

    smooth_limit = max(
        2.0, symmetric_limit(smoothed_standardized, 99.0)
    )
    smooth_image = axes[2].imshow(
        smoothed_standardized,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=[
            frequency_mhz[0],
            frequency_mhz[-1],
            binned_elapsed_minutes[0],
            binned_elapsed_minutes[-1],
        ],
        cmap="RdBu_r",
        vmin=-smooth_limit,
        vmax=smooth_limit,
    )
    axes[2].set_ylabel("time (min)")
    axes[2].set_xlabel("topocentric frequency (MHz)")
    axes[2].set_title(
        f"{args.smooth_fwhm_mhz:.3f} MHz frequency-smoothed response"
    )
    figure.colorbar(
        smooth_image,
        ax=axes[2],
        pad=0.01,
        label="smoothed response / robust RMS",
    )
    figure.suptitle(
        f"{args.sample_id}: pre-Doppler S/{args.field}, "
        f"{args.polar}"
    )
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=170)
    plt.close(figure)

    print(f"sample_id={args.sample_id}")
    print(f"selected_spectra={len(values)}")
    print(f"time_bins={len(binned_values)}")
    print(f"time_bin_seconds={samples_per_bin * sample_seconds:.6f}")
    print(f"existing_rfi_fraction={np.mean(existing_rfi):.9f}")
    print(f"existing_excluded_fraction={np.mean(existing_excluded):.9f}")
    print(f"saved_npz={args.output_npz}")
    print(f"saved_figure={args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
