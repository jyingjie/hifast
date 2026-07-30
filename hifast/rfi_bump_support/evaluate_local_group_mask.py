#!/usr/bin/env python3
"""Fit and evaluate masks for one beam-local dynamically tracked comb group."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import shift as shift_row

from build_candidate_dynamic_masks import gaussian_coverage, group_mask
from evaluate_candidate_masks import (
    profile_leakage,
    weighted_aligned_mean,
)
from fit_comb_track_profiles import (
    fit_profile,
    gaussian_with_linear_background,
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--waterfall", type=Path, required=True)
    parser.add_argument("--primary-tracks", type=Path, required=True)
    parser.add_argument("--comparison-tracks", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--group-label", default="主组")
    parser.add_argument("--output-masks", type=Path, required=True)
    parser.add_argument("--profile-csv", type=Path, required=True)
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--tooth-metrics-csv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--fit-half-width-mhz", type=float, default=0.12)
    parser.add_argument("--center-bound-mhz", type=float, default=0.06)
    parser.add_argument("--core-half-width-fwhm", type=float, default=0.5)
    parser.add_argument(
        "--conservative-half-width-fwhm", type=float, default=0.75
    )
    parser.add_argument(
        "--conservative-margin-channels", type=float, default=1.0
    )
    parser.add_argument(
        "--maximum-track-disagreement-channels",
        type=float,
        default=4.0,
    )
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def binned_flag_to_raw(
    binned_flag: np.ndarray,
    starts: np.ndarray,
    stops: np.ndarray,
    raw_count: int,
) -> np.ndarray:
    result = np.zeros(raw_count, dtype=bool)
    for flag, start, stop in zip(binned_flag, starts, stops):
        result[start:stop] = flag
    return result


def main() -> int:
    arguments = parse_arguments()
    with np.load(arguments.waterfall, allow_pickle=False) as product:
        frequency_mhz = product["frequency_mhz"].astype(np.float64)
        raw_residual = product["raw_residual"].astype(np.float64)
        binned_residual = product["binned_residual"].astype(np.float64)
        selected_mjd = product["selected_mjd"].astype(np.float64)
        binned_mjd = product["binned_mjd"].astype(np.float64)
        bin_starts = product["bin_start_offsets"].astype(np.int64)
        bin_stops = product["bin_stop_offsets"].astype(np.int64)
        selected_original_indices = product[
            "selected_original_indices"
        ].astype(np.int64)
        source_hdf5 = str(product["source_hdf5"])
        polarization = str(product["polarization"])
        data_unit = str(product["data_unit"])
    with np.load(arguments.primary_tracks, allow_pickle=False) as product:
        guides_mhz = product["guide_centers_mhz"].astype(np.float64)
        primary_mjd = product["binned_mjd"].astype(np.float64)
        primary_shift_mhz = product[
            "independent_common_shift_mhz"
        ].astype(np.float64)
        primary_track_mhz = product[
            "independent_common_track_mhz"
        ].astype(np.float64)
        primary_reliable = product[
            "independent_reliable"
        ].astype(bool)
    with np.load(arguments.comparison_tracks, allow_pickle=False) as product:
        comparison_mjd = product["binned_mjd"].astype(np.float64)
        comparison_track_mhz = product[
            "independent_common_track_mhz"
        ].astype(np.float64)

    if len(primary_mjd) != len(binned_mjd) or not np.allclose(
        primary_mjd, binned_mjd, rtol=0, atol=1e-10
    ):
        raise ValueError(
            "Primary tracks do not use the same time bins as the Waterfall"
        )
    if primary_track_mhz.shape[1] != comparison_track_mhz.shape[1]:
        raise ValueError(
            "Primary and comparison tracks have different tooth counts"
        )

    delta_frequency_mhz = float(np.median(np.diff(frequency_mhz)))
    primary_shift_channels = primary_shift_mhz / delta_frequency_mhz
    aligned = np.empty_like(binned_residual)
    for time_index, (row, shift_channels) in enumerate(
        zip(binned_residual, primary_shift_channels)
    ):
        aligned[time_index] = shift_row(
            row,
            -shift_channels,
            order=1,
            mode="nearest",
            prefilter=False,
        )
    stacked_profile = np.nanmean(aligned, axis=0)

    profile_rows: list[dict[str, object]] = []
    fitted_parameters = []
    fitted_centers_mhz = []
    fitted_fwhm_mhz = []
    fitted_models = []
    for tooth_number, guide_mhz in enumerate(guides_mhz, start=1):
        fitted, covariance, noise = fit_profile(
            frequency_mhz,
            stacked_profile,
            float(guide_mhz),
            arguments.fit_half_width_mhz,
            arguments.center_bound_mhz,
        )
        center_mhz = float(fitted[3])
        sigma_mhz = float(fitted[4])
        fwhm_mhz = 2.354820045 * sigma_mhz
        amplitude = float(fitted[2])
        center_uncertainty_mhz = float(
            np.sqrt(max(0.0, covariance[3, 3]))
        )
        profile_rows.append(
            {
                "tooth_number": tooth_number,
                "guide_center_mhz": f"{guide_mhz:.9f}",
                "fitted_center_mhz": f"{center_mhz:.9f}",
                "center_offset_khz": (
                    f"{(center_mhz - guide_mhz) * 1000.0:.6f}"
                ),
                "center_uncertainty_khz": (
                    f"{center_uncertainty_mhz * 1000.0:.6f}"
                ),
                "fwhm_khz": f"{fwhm_mhz * 1000.0:.6f}",
                "amplitude": f"{amplitude:.12g}",
                "residual_noise": f"{noise:.12g}",
                "amplitude_snr": (
                    f"{amplitude / noise:.6f}" if noise > 0 else "nan"
                ),
                "data_unit": data_unit,
            }
        )
        fitted_parameters.append(fitted)
        fitted_centers_mhz.append(center_mhz)
        fitted_fwhm_mhz.append(fwhm_mhz)
        selected = (
            np.abs(frequency_mhz - float(guide_mhz))
            <= arguments.fit_half_width_mhz
        )
        fitted_models.append(
            (
                frequency_mhz[selected],
                gaussian_with_linear_background(
                    frequency_mhz[selected], *fitted
                ),
            )
        )

    fitted_parameters_array = np.asarray(fitted_parameters)
    fitted_centers_array = np.asarray(fitted_centers_mhz)
    fitted_fwhm_array = np.asarray(fitted_fwhm_mhz)
    write_csv(arguments.profile_csv, profile_rows)

    comparison_at_primary = np.column_stack(
        [
            np.interp(
                primary_mjd,
                comparison_mjd,
                comparison_track_mhz[:, tooth_index],
            )
            for tooth_index in range(comparison_track_mhz.shape[1])
        ]
    )
    track_difference_channels = (
        np.abs(primary_track_mhz - comparison_at_primary)
        / abs(delta_frequency_mhz)
    )
    agreement_binned = (
        np.max(track_difference_channels, axis=1)
        <= arguments.maximum_track_disagreement_channels
    )
    accepted_binned = primary_reliable & agreement_binned
    accepted_raw = binned_flag_to_raw(
        accepted_binned, bin_starts, bin_stops, len(selected_mjd)
    )

    raw_shift_mhz = np.interp(
        selected_mjd, primary_mjd, primary_shift_mhz
    )
    core_half_widths = (
        arguments.core_half_width_fwhm * fitted_fwhm_array
    )
    conservative_half_widths = (
        arguments.conservative_half_width_fwhm * fitted_fwhm_array
        + arguments.conservative_margin_channels
        * abs(delta_frequency_mhz)
    )
    core_all = group_mask(
        frequency_mhz,
        raw_shift_mhz,
        fitted_centers_array,
        core_half_widths,
    )
    conservative_all = group_mask(
        frequency_mhz,
        raw_shift_mhz,
        fitted_centers_array,
        conservative_half_widths,
    )
    masks = {
        "core_all_times": core_all,
        "conservative_all_times": conservative_all,
        "core_confidence_gated": core_all & accepted_raw[:, None],
        "conservative_confidence_gated": (
            conservative_all & accepted_raw[:, None]
        ),
    }

    raw_shift_channels = raw_shift_mhz / delta_frequency_mhz
    original_valid = np.isfinite(raw_residual)
    original_profile, original_weight = weighted_aligned_mean(
        raw_residual, original_valid, raw_shift_channels
    )
    clean_frequency = np.ones(len(frequency_mhz), dtype=bool)
    for center_mhz in fitted_centers_array:
        clean_frequency &= np.abs(
            frequency_mhz - center_mhz
        ) > arguments.fit_half_width_mhz

    metric_rows: list[dict[str, object]] = []
    tooth_rows: list[dict[str, object]] = []
    masked_profiles = {}
    masked_weights = {}
    core_coverage = gaussian_coverage(
        fitted_fwhm_array, core_half_widths
    )
    conservative_coverage = gaussian_coverage(
        fitted_fwhm_array, conservative_half_widths
    )
    for candidate, mask in masks.items():
        masked_valid = original_valid & ~mask
        masked_profile, masked_weight = weighted_aligned_mean(
            raw_residual, masked_valid, raw_shift_channels
        )
        leakage, per_tooth = profile_leakage(
            frequency_mhz,
            original_profile,
            original_weight,
            masked_profile,
            masked_weight,
            fitted_parameters_array,
        )
        masked_profiles[candidate] = masked_profile
        masked_weights[candidate] = masked_weight
        row_fraction = np.mean(mask, axis=1)
        channel_time_loss = np.mean(mask, axis=0)
        clean = (
            clean_frequency
            & np.isfinite(original_profile)
            & np.isfinite(masked_profile)
        )
        clean_difference = masked_profile[clean] - original_profile[clean]
        geometric_coverage = (
            core_coverage
            if candidate.startswith("core")
            else conservative_coverage
        )
        metric_rows.append(
            {
                "candidate": candidate,
                "masked_fraction": f"{np.mean(mask):.9f}",
                "active_time_fraction": (
                    f"{np.mean(np.any(mask, axis=1)):.9f}"
                ),
                "row_mask_fraction_median": (
                    f"{np.median(row_fraction):.9f}"
                ),
                "channel_time_loss_p95": (
                    f"{np.percentile(channel_time_loss, 95):.9f}"
                ),
                "channel_time_loss_max": (
                    f"{np.max(channel_time_loss):.9f}"
                ),
                "geometric_coverage_mean": (
                    f"{np.mean(geometric_coverage):.9f}"
                ),
                "estimated_leakage_fraction": f"{leakage:.9f}",
                "estimated_coverage_fraction": f"{1.0 - leakage:.9f}",
                "clean_profile_difference_median": (
                    f"{np.median(clean_difference):.12g}"
                ),
                "clean_profile_difference_rms": (
                    f"{np.sqrt(np.mean(clean_difference**2)):.12g}"
                ),
                "outside_mask_values_modified": "false",
            }
        )
        for tooth_number, tooth_leakage in enumerate(
            per_tooth, start=1
        ):
            tooth_rows.append(
                {
                    "candidate": candidate,
                    "tooth_number": tooth_number,
                    "estimated_leakage_fraction": (
                        f"{tooth_leakage:.9f}"
                    ),
                }
            )
    write_csv(arguments.metrics_csv, metric_rows)
    write_csv(arguments.tooth_metrics_csv, tooth_rows)

    arguments.output_masks.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        arguments.output_masks,
        sample_id=arguments.sample_id,
        source_hdf5=source_hdf5,
        source_waterfall=str(arguments.waterfall),
        source_primary_tracks=str(arguments.primary_tracks),
        source_comparison_tracks=str(arguments.comparison_tracks),
        polarization=polarization,
        frequency_mhz=frequency_mhz,
        selected_original_indices=selected_original_indices,
        selected_mjd=selected_mjd,
        raw_common_shift_mhz=raw_shift_mhz,
        guide_centers_mhz=guides_mhz,
        fitted_parameters=fitted_parameters_array,
        fitted_centers_mhz=fitted_centers_array,
        fitted_fwhm_mhz=fitted_fwhm_array,
        core_half_widths_mhz=core_half_widths,
        conservative_half_widths_mhz=conservative_half_widths,
        track_difference_channels=track_difference_channels,
        track_agreement_binned=agreement_binned,
        primary_reliable_binned=primary_reliable,
        accepted_binned=accepted_binned,
        accepted_raw=accepted_raw,
        mask_core_all_times=masks["core_all_times"],
        mask_conservative_all_times=masks[
            "conservative_all_times"
        ],
        mask_core_confidence_gated=masks[
            "core_confidence_gated"
        ],
        mask_conservative_confidence_gated=masks[
            "conservative_confidence_gated"
        ],
        formal_science_mask=False,
    )

    scale = 1000.0 if data_unit in {"Jy/beam", "K"} else 1.0
    display_unit = (
        "mJy/beam"
        if data_unit == "Jy/beam"
        else "mK"
        if data_unit == "K"
        else data_unit
    )
    figure, axes = plt.subplots(
        3, 1, figsize=(15, 10), constrained_layout=True
    )
    axes[0].plot(
        frequency_mhz,
        stacked_profile * scale,
        color="0.45",
        linewidth=0.9,
        label="de-drifted mean",
    )
    for tooth_number, (x, model) in enumerate(
        fitted_models, start=1
    ):
        axes[0].plot(
            x,
            model * scale,
            linewidth=1.2,
            label=f"tooth {tooth_number}",
        )
    axes[0].set(
        ylabel=f"residual ({display_unit})",
        title=f"{arguments.sample_id}: local profile fits",
    )
    axes[0].legend(ncol=len(fitted_models) + 1, fontsize=8)
    axes[0].grid(alpha=0.2)

    axes[1].plot(
        frequency_mhz,
        original_profile * scale,
        color="0.4",
        linewidth=0.8,
        label="original",
    )
    for candidate, color in (
        ("core_all_times", "tab:blue"),
        ("conservative_all_times", "tab:red"),
        ("conservative_confidence_gated", "tab:green"),
    ):
        axes[1].plot(
            frequency_mhz,
            masked_profiles[candidate] * scale,
            color=color,
            linewidth=0.8,
            label=candidate,
        )
    axes[1].set(
        ylabel=f"residual ({display_unit})",
        title="Available de-drifted profile after masking",
    )
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.2)

    for candidate, color in (
        ("core_all_times", "tab:blue"),
        ("conservative_all_times", "tab:red"),
        ("conservative_confidence_gated", "tab:green"),
    ):
        axes[2].plot(
            frequency_mhz,
            np.mean(masks[candidate], axis=0) * 100.0,
            color=color,
            linewidth=1.0,
            label=candidate,
        )
    axes[2].set(
        xlabel="de-drifted topocentric frequency (MHz)",
        ylabel="time samples masked (%)",
        title="Per-channel time loss",
    )
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.2)
    arguments.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.figure, dpi=180)
    plt.close(figure)

    lines = [
        f"# {arguments.sample_id} 局部{arguments.group_label}候选 mask 评价",
        "",
        f"本报告只评价当前局部{arguments.group_label}。",
        "",
        "## 轨迹和宽度",
        "",
        f"- 15/30 秒轨迹全部 {len(guides_mhz)} 齿相差不超过 "
        f"`{arguments.maximum_track_disagreement_channels:g} channel` "
        f"的时间比例：`{np.mean(agreement_binned) * 100:.2f}%`；",
        f"- 30 秒轨迹自身可靠时间比例："
        f"`{np.mean(primary_reliable) * 100:.2f}%`；",
        f"- 同时满足轨迹可靠和跨时间合并一致的比例："
        f"`{np.mean(accepted_binned) * 100:.2f}%`；",
        f"- {len(guides_mhz)} 齿 FWHM："
        + "、".join(
            f"`{value * 1000.0:.2f} kHz`"
            for value in fitted_fwhm_array
        )
        + "。",
        "",
        "## 候选结果",
        "",
        "| 候选 | 遮挡比例 | 估计残留 | 最大单通道时间损失 |",
        "|---|---:|---:|---:|",
    ]
    for row in metric_rows:
        lines.append(
            f"| {row['candidate']} | "
            f"{float(row['masked_fraction']) * 100:.2f}% | "
            f"{float(row['estimated_leakage_fraction']) * 100:.2f}% | "
            f"{float(row['channel_time_loss_max']) * 100:.2f}% |"
        )
    lines.extend(
        [
            "",
            "所有 mask 都是独立测试数组；原始 HDF5 未修改，"
            "也没有生成 cube。",
            "",
        ]
    )
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(
        "\n".join(lines), encoding="utf-8"
    )

    print(f"sample_id={arguments.sample_id}")
    print(
        "track_agreement_fraction="
        f"{np.mean(agreement_binned):.9f}"
    )
    print(
        "accepted_time_fraction="
        f"{np.mean(accepted_binned):.9f}"
    )
    print(
        "fitted_fwhm_khz="
        + ",".join(
            f"{value * 1000.0:.6f}" for value in fitted_fwhm_array
        )
    )
    for row in metric_rows:
        print(
            f"{row['candidate']}: "
            f"mask={row['masked_fraction']},"
            f"leakage={row['estimated_leakage_fraction']},"
            f"channel_loss_max={row['channel_time_loss_max']}"
        )
    print(f"saved_masks={arguments.output_masks}")
    print(f"saved_report={arguments.report}")
    print(f"saved_figure={arguments.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
