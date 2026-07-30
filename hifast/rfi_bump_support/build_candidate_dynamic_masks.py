#!/usr/bin/env python3
"""Build reversible candidate masks from two shared-shift comb groups."""

from __future__ import annotations

import argparse
import csv
from math import erf, sqrt
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def group_mask(
    frequency_mhz: np.ndarray,
    raw_shift_mhz: np.ndarray,
    centers_mhz: np.ndarray,
    half_widths_mhz: np.ndarray,
) -> np.ndarray:
    mask = np.zeros(
        (len(raw_shift_mhz), len(frequency_mhz)), dtype=bool
    )
    for center_mhz, half_width_mhz in zip(
        centers_mhz, half_widths_mhz
    ):
        dynamic_center = center_mhz + raw_shift_mhz
        mask |= (
            np.abs(
                frequency_mhz[None, :] - dynamic_center[:, None]
            )
            <= half_width_mhz
        )
    return mask


def gaussian_coverage(
    fwhm_mhz: np.ndarray,
    half_width_mhz: np.ndarray,
) -> np.ndarray:
    sigma_mhz = fwhm_mhz / 2.354820045
    return np.asarray(
        [
            erf(float(width / sigma) / sqrt(2))
            for width, sigma in zip(half_width_mhz, sigma_mhz)
        ]
    )


def write_metrics(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def binned_mask(
    raw_mask: np.ndarray,
    starts: np.ndarray,
    stops: np.ndarray,
) -> np.ndarray:
    output = np.zeros((len(starts), raw_mask.shape[1]), dtype=bool)
    for index, (start, stop) in enumerate(zip(starts, stops)):
        output[index] = np.any(raw_mask[start:stop], axis=0)
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--waterfall", type=Path, required=True)
    parser.add_argument("--group-a-tracks", type=Path, required=True)
    parser.add_argument("--profile-fits", type=Path, required=True)
    parser.add_argument("--output-npz", type=Path, required=True)
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--core-half-width-fwhm", type=float, default=0.5)
    parser.add_argument(
        "--conservative-half-width-fwhm", type=float, default=0.75
    )
    parser.add_argument(
        "--conservative-margin-channels", type=float, default=1.0
    )
    parser.add_argument(
        "--group-b-support-threshold", type=float, default=0.5
    )
    parser.add_argument(
        "--group-b-minimum-supported-teeth", type=int, default=2
    )
    args = parser.parse_args()

    with np.load(args.waterfall, allow_pickle=False) as product:
        frequency_mhz = product["frequency_mhz"].astype(np.float64)
        raw_elapsed_minutes = product["elapsed_minutes"].astype(np.float64)
        binned_elapsed_minutes = product[
            "binned_elapsed_minutes"
        ].astype(np.float64)
        waterfall = product["smoothed_standardized"].astype(np.float64)
        selected_original_indices = product[
            "selected_original_indices"
        ].astype(np.int64)
        selected_mjd = product["selected_mjd"].astype(np.float64)
        bin_start_offsets = product["bin_start_offsets"].astype(np.int64)
        bin_stop_offsets = product["bin_stop_offsets"].astype(np.int64)
        source_hdf5 = str(product["source_hdf5"])
        polarization = str(product["polarization"])
    with np.load(args.group_a_tracks, allow_pickle=False) as product:
        common_shift_mhz = product["common_shift_mhz"].astype(np.float64)
    with np.load(args.profile_fits, allow_pickle=False) as product:
        fitted_groups = product["fitted_groups"].astype(str)
        fitted_centers_mhz = product[
            "fitted_centers_mhz"
        ].astype(np.float64)
        fitted_fwhm_mhz = product["fitted_fwhm_mhz"].astype(np.float64)

    delta_frequency_mhz = abs(float(np.median(np.diff(frequency_mhz))))
    raw_shift_mhz = np.interp(
        raw_elapsed_minutes,
        binned_elapsed_minutes,
        common_shift_mhz,
    )
    is_group_a = fitted_groups == "A_stronger"
    is_group_b = fitted_groups == "B_weaker"
    centers_a = fitted_centers_mhz[is_group_a]
    centers_b = fitted_centers_mhz[is_group_b]
    fwhm_a = fitted_fwhm_mhz[is_group_a]
    fwhm_b = fitted_fwhm_mhz[is_group_b]

    core_half_a = args.core_half_width_fwhm * fwhm_a
    core_half_b = args.core_half_width_fwhm * fwhm_b
    conservative_margin_mhz = (
        args.conservative_margin_channels * delta_frequency_mhz
    )
    conservative_half_a = (
        args.conservative_half_width_fwhm * fwhm_a
        + conservative_margin_mhz
    )
    conservative_half_b = (
        args.conservative_half_width_fwhm * fwhm_b
        + conservative_margin_mhz
    )

    core_a = group_mask(
        frequency_mhz, raw_shift_mhz, centers_a, core_half_a
    )
    core_b = group_mask(
        frequency_mhz, raw_shift_mhz, centers_b, core_half_b
    )
    conservative_a = group_mask(
        frequency_mhz,
        raw_shift_mhz,
        centers_a,
        conservative_half_a,
    )
    conservative_b = group_mask(
        frequency_mhz,
        raw_shift_mhz,
        centers_b,
        conservative_half_b,
    )

    shared_b_tracks = (
        centers_b[None, :] + common_shift_mhz[:, None]
    )
    shared_b_channels = np.abs(
        frequency_mhz[None, None, :]
        - shared_b_tracks[:, :, None]
    ).argmin(axis=2)
    time_indices = np.arange(len(waterfall))[:, None]
    shared_b_responses = waterfall[time_indices, shared_b_channels]
    b_supported_teeth = np.count_nonzero(
        shared_b_responses >= args.group_b_support_threshold,
        axis=1,
    )
    b_active_binned = (
        b_supported_teeth >= args.group_b_minimum_supported_teeth
    )
    b_active_raw = np.zeros(len(raw_elapsed_minutes), dtype=bool)
    for active, start, stop in zip(
        b_active_binned, bin_start_offsets, bin_stop_offsets
    ):
        b_active_raw[start:stop] = active

    candidates = {
        "core_all_groups": core_a | core_b,
        "conservative_all_groups": conservative_a | conservative_b,
        "conservative_a_only": conservative_a,
        "conservative_a_plus_gated_b": (
            conservative_a | (conservative_b & b_active_raw[:, None])
        ),
    }
    coverage_core_a = gaussian_coverage(fwhm_a, core_half_a)
    coverage_core_b = gaussian_coverage(fwhm_b, core_half_b)
    coverage_conservative_a = gaussian_coverage(
        fwhm_a, conservative_half_a
    )
    coverage_conservative_b = gaussian_coverage(
        fwhm_b, conservative_half_b
    )
    b_active_fraction = float(np.mean(b_active_raw))
    coverage_by_candidate = {
        "core_all_groups": (
            float(np.mean(coverage_core_a)),
            float(np.mean(coverage_core_b)),
        ),
        "conservative_all_groups": (
            float(np.mean(coverage_conservative_a)),
            float(np.mean(coverage_conservative_b)),
        ),
        "conservative_a_only": (
            float(np.mean(coverage_conservative_a)),
            0.0,
        ),
        "conservative_a_plus_gated_b": (
            float(np.mean(coverage_conservative_a)),
            float(np.mean(coverage_conservative_b)) * b_active_fraction,
        ),
    }

    rows = []
    for candidate_name, mask in candidates.items():
        coverage_a, coverage_b = coverage_by_candidate[candidate_name]
        rows.append(
            {
                "candidate": candidate_name,
                "masked_fraction_test_band": f"{np.mean(mask):.9f}",
                "masked_fraction_a": (
                    f"{np.mean(mask & conservative_a):.9f}"
                    if "conservative" in candidate_name
                    else f"{np.mean(mask & core_a):.9f}"
                ),
                "masked_fraction_b": (
                    f"{np.mean(mask & conservative_b):.9f}"
                    if "conservative" in candidate_name
                    else f"{np.mean(mask & core_b):.9f}"
                ),
                "model_energy_coverage_a": f"{coverage_a:.9f}",
                "model_energy_coverage_b": f"{coverage_b:.9f}",
                "group_b_active_time_fraction": (
                    f"{b_active_fraction:.9f}"
                ),
            }
        )
    write_metrics(args.metrics_csv, rows)

    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        source_hdf5=source_hdf5,
        source_waterfall=str(args.waterfall),
        source_tracks=str(args.group_a_tracks),
        source_profile_fits=str(args.profile_fits),
        polarization=polarization,
        frequency_mhz=frequency_mhz,
        selected_original_indices=selected_original_indices,
        selected_mjd=selected_mjd,
        raw_elapsed_minutes=raw_elapsed_minutes,
        raw_common_shift_mhz=raw_shift_mhz,
        fitted_groups=fitted_groups,
        fitted_centers_mhz=fitted_centers_mhz,
        fitted_fwhm_mhz=fitted_fwhm_mhz,
        core_half_widths_mhz=np.concatenate([core_half_a, core_half_b]),
        conservative_half_widths_mhz=np.concatenate(
            [conservative_half_a, conservative_half_b]
        ),
        group_b_active_binned=b_active_binned,
        group_b_active_raw=b_active_raw,
        mask_core_all_groups=candidates["core_all_groups"],
        mask_conservative_all_groups=candidates[
            "conservative_all_groups"
        ],
        mask_conservative_a_only=candidates["conservative_a_only"],
        mask_conservative_a_plus_gated_b=candidates[
            "conservative_a_plus_gated_b"
        ],
        formal_science_mask=False,
    )

    binned_core = binned_mask(
        candidates["core_all_groups"],
        bin_start_offsets,
        bin_stop_offsets,
    )
    binned_gated = binned_mask(
        candidates["conservative_a_plus_gated_b"],
        bin_start_offsets,
        bin_stop_offsets,
    )
    figure, axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True)
    limit = max(2.0, float(np.nanpercentile(np.abs(waterfall), 99)))
    image = axes[0].imshow(
        waterfall,
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
        vmin=-limit,
        vmax=limit,
    )
    axes[0].set_ylabel("time (min)")
    axes[0].set_title("Detection Waterfall")
    figure.colorbar(image, ax=axes[0], pad=0.01, label="response")
    for axis, mask, title in (
        (axes[1], binned_core, "Core mask: all A and B teeth"),
        (
            axes[2],
            binned_gated,
            "Conservative A plus confidence-gated B",
        ),
    ):
        axis.imshow(
            mask,
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            extent=[
                frequency_mhz[0],
                frequency_mhz[-1],
                binned_elapsed_minutes[0],
                binned_elapsed_minutes[-1],
            ],
            cmap="gray_r",
            vmin=0,
            vmax=1,
        )
        axis.set_ylabel("time (min)")
        axis.set_title(title)
    axes[2].set_xlabel("topocentric frequency (MHz)")
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=170)
    plt.close(figure)

    for row in rows:
        print(
            f"{row['candidate']}: "
            f"masked_fraction={row['masked_fraction_test_band']},"
            f"coverage_A={row['model_energy_coverage_a']},"
            f"coverage_B={row['model_energy_coverage_b']}"
        )
    print(f"saved_npz={args.output_npz}")
    print(f"saved_metrics={args.metrics_csv}")
    print(f"saved_figure={args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
