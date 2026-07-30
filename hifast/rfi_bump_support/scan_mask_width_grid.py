#!/usr/bin/env python3
"""Scan conservative mask widths across track-consistent beams."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from build_candidate_dynamic_masks import group_mask
from evaluate_candidate_masks import profile_leakage, weighted_aligned_mean


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--track-consistency-csv", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--scan", type=int, default=43)
    parser.add_argument("--polar", choices=["XX", "YY"], default="YY")
    parser.add_argument("--beams", nargs="+")
    parser.add_argument(
        "--waterfall-variant", default="pre_fc_30s"
    )
    parser.add_argument("--fit-label", required=True)
    parser.add_argument(
        "--half-width-fwhm", type=float, nargs="+", required=True
    )
    parser.add_argument(
        "--margin-channels", type=float, nargs="+", required=True
    )
    parser.add_argument(
        "--minimum-amplitude-snr", type=float, nargs="+", default=[0.0]
    )
    parser.add_argument(
        "--maximum-fwhm-khz", type=float, default=float("inf")
    )
    parser.add_argument(
        "--minimum-masked-tooth-count", type=int, default=3
    )
    parser.add_argument(
        "--maximum-leakage-fraction", type=float, default=0.10
    )
    parser.add_argument(
        "--maximum-tooth-leakage-fraction", type=float, default=0.15
    )
    parser.add_argument(
        "--maximum-mask-fraction", type=float, default=0.15
    )
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-selection-csv", type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    arguments = parse_arguments()
    beams = [
        row["beam"]
        for row in read_csv(arguments.track_consistency_csv)
        if row["track_consistency_accepted"].lower() == "true"
        and (
            arguments.beams is None
            or row["beam"] in arguments.beams
        )
    ]
    rows = []
    beam_diagnostics = {}
    for beam in beams:
        waterfall_path = arguments.tables_dir / (
            f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
            f"{arguments.waterfall_variant}_waterfall.npz"
        )
        fit_path = arguments.tables_dir / (
            f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
            f"{arguments.fit_label}_candidate_masks.npz"
        )
        profile_path = arguments.tables_dir / (
            f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
            f"{arguments.fit_label}_profile_fits.csv"
        )
        with np.load(waterfall_path, allow_pickle=False) as product:
            frequency_mhz = product["frequency_mhz"].astype(np.float64)
            raw_residual = product["raw_residual"].astype(np.float64)
        with np.load(fit_path, allow_pickle=False) as product:
            raw_shift_mhz = product[
                "raw_common_shift_mhz"
            ].astype(np.float64)
            fitted_parameters = product[
                "fitted_parameters"
            ].astype(np.float64)
            fitted_centers = product[
                "fitted_centers_mhz"
            ].astype(np.float64)
            fitted_fwhm = product[
                "fitted_fwhm_mhz"
            ].astype(np.float64)
        amplitude_snr = np.asarray(
            [
                float(row["amplitude_snr"])
                for row in read_csv(profile_path)
            ],
            dtype=np.float64,
        )
        beam_diagnostics[beam] = {
            "total_tooth_count": len(amplitude_snr),
            "amplitude_snr": amplitude_snr,
            "fitted_fwhm_khz": fitted_fwhm * 1000.0,
        }
        delta_frequency_mhz = float(np.median(np.diff(frequency_mhz)))
        shift_channels = raw_shift_mhz / delta_frequency_mhz
        original_valid = np.isfinite(raw_residual)
        original_profile, original_weight = weighted_aligned_mean(
            raw_residual, original_valid, shift_channels
        )
        for half_width_fwhm in arguments.half_width_fwhm:
            for margin_channels in arguments.margin_channels:
                for minimum_amplitude_snr in (
                    arguments.minimum_amplitude_snr
                ):
                    selected_teeth = (
                        amplitude_snr >= minimum_amplitude_snr
                    ) & (
                        fitted_fwhm * 1000.0
                        <= arguments.maximum_fwhm_khz
                    )
                    if (
                        np.count_nonzero(selected_teeth)
                        < arguments.minimum_masked_tooth_count
                    ):
                        continue
                    half_widths = (
                        half_width_fwhm * fitted_fwhm[selected_teeth]
                        + margin_channels * abs(delta_frequency_mhz)
                    )
                    mask = group_mask(
                        frequency_mhz,
                        raw_shift_mhz,
                        fitted_centers[selected_teeth],
                        half_widths,
                    )
                    masked_profile, masked_weight = weighted_aligned_mean(
                        raw_residual,
                        original_valid & ~mask,
                        shift_channels,
                    )
                    leakage, per_tooth = profile_leakage(
                        frequency_mhz,
                        original_profile,
                        original_weight,
                        masked_profile,
                        masked_weight,
                        fitted_parameters[selected_teeth],
                    )
                    per_tooth_array = np.asarray(
                        per_tooth, dtype=np.float64
                    )
                    all_tooth_metrics_finite = bool(
                        np.all(np.isfinite(per_tooth_array))
                    )
                    maximum_tooth_leakage = (
                        float(np.max(per_tooth_array))
                        if all_tooth_metrics_finite
                        else np.nan
                    )
                    channel_time_loss = np.mean(mask, axis=0)
                    rows.append(
                        {
                            "scan": arguments.scan,
                            "beam": beam,
                            "polar": arguments.polar,
                            "half_width_fwhm": half_width_fwhm,
                            "margin_channels": margin_channels,
                            "minimum_amplitude_snr": minimum_amplitude_snr,
                            "maximum_fwhm_khz": (
                                arguments.maximum_fwhm_khz
                            ),
                            "masked_tooth_count": int(
                                np.count_nonzero(selected_teeth)
                            ),
                            "total_tooth_count": len(selected_teeth),
                            "masked_fraction": float(np.mean(mask)),
                            "estimated_leakage_fraction": leakage,
                            "maximum_tooth_leakage_fraction": (
                                maximum_tooth_leakage
                            ),
                            "all_tooth_leakage_metrics_finite": (
                                all_tooth_metrics_finite
                            ),
                            "channel_time_loss_p95": float(
                                np.percentile(channel_time_loss, 95)
                            ),
                            "channel_time_loss_max": float(
                                np.max(channel_time_loss)
                            ),
                        }
                    )
    grid_fieldnames = [
        "scan",
        "beam",
        "polar",
        "half_width_fwhm",
        "margin_channels",
        "minimum_amplitude_snr",
        "maximum_fwhm_khz",
        "masked_tooth_count",
        "total_tooth_count",
        "masked_fraction",
        "estimated_leakage_fraction",
        "maximum_tooth_leakage_fraction",
        "all_tooth_leakage_metrics_finite",
        "channel_time_loss_p95",
        "channel_time_loss_max",
    ]
    arguments.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output_csv.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=grid_fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    if arguments.output_selection_csv is not None:
        selected_rows = []
        for beam in beams:
            eligible = [
                row
                for row in rows
                if row["beam"] == beam
                and row["all_tooth_leakage_metrics_finite"]
                and row["masked_fraction"]
                <= arguments.maximum_mask_fraction
                and row["estimated_leakage_fraction"]
                <= arguments.maximum_leakage_fraction
                and row["maximum_tooth_leakage_fraction"]
                <= arguments.maximum_tooth_leakage_fraction
            ]
            eligible.sort(
                key=lambda row: (
                    row["masked_fraction"],
                    row["estimated_leakage_fraction"],
                    row["half_width_fwhm"],
                    row["margin_channels"],
                )
            )
            if eligible:
                selected = {
                    **eligible[0],
                    "adaptive_mask_accepted": True,
                    "selection_rejection_reason": "",
                }
            else:
                beam_rows = [
                    row for row in rows if row["beam"] == beam
                ]
                if beam_rows:
                    selected = {
                        **min(
                            beam_rows,
                            key=lambda row: (
                                max(
                                    0.0,
                                    row["masked_fraction"]
                                    - arguments.maximum_mask_fraction,
                                )
                                + max(
                                    0.0,
                                    row["estimated_leakage_fraction"]
                                    - arguments.maximum_leakage_fraction,
                                ),
                                row["masked_fraction"],
                            ),
                        ),
                        "adaptive_mask_accepted": False,
                        "selection_rejection_reason": (
                            "mask_or_leakage_limits_not_met"
                        ),
                    }
                else:
                    diagnostic = beam_diagnostics[beam]
                    eligible_teeth = (
                        diagnostic["amplitude_snr"]
                        >= min(arguments.minimum_amplitude_snr)
                    ) & (
                        diagnostic["fitted_fwhm_khz"]
                        <= arguments.maximum_fwhm_khz
                    )
                    selected = {
                        "scan": arguments.scan,
                        "beam": beam,
                        "polar": arguments.polar,
                        "half_width_fwhm": "",
                        "margin_channels": min(
                            arguments.margin_channels
                        ),
                        "minimum_amplitude_snr": min(
                            arguments.minimum_amplitude_snr
                        ),
                        "maximum_fwhm_khz": (
                            arguments.maximum_fwhm_khz
                        ),
                        "masked_tooth_count": int(
                            np.count_nonzero(eligible_teeth)
                        ),
                        "total_tooth_count": diagnostic[
                            "total_tooth_count"
                        ],
                        "masked_fraction": "",
                        "estimated_leakage_fraction": "",
                        "maximum_tooth_leakage_fraction": "",
                        "all_tooth_leakage_metrics_finite": False,
                        "channel_time_loss_p95": "",
                        "channel_time_loss_max": "",
                        "adaptive_mask_accepted": False,
                        "selection_rejection_reason": (
                            "fewer_than_minimum_shape_accepted_teeth"
                        ),
                    }
            selected_rows.append(selected)
        arguments.output_selection_csv.parent.mkdir(
            parents=True, exist_ok=True
        )
        with arguments.output_selection_csv.open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(selected_rows[0])
            )
            writer.writeheader()
            writer.writerows(selected_rows)
        accepted = sum(
            row["adaptive_mask_accepted"] for row in selected_rows
        )
        print(
            f"adaptive_pass={accepted}/{len(selected_rows)}, "
            f"saved_selection_csv={arguments.output_selection_csv}"
        )
    for half_width_fwhm in arguments.half_width_fwhm:
        for margin_channels in arguments.margin_channels:
            for minimum_amplitude_snr in (
                arguments.minimum_amplitude_snr
            ):
                selected = [
                    row
                    for row in rows
                    if row["half_width_fwhm"] == half_width_fwhm
                    and row["margin_channels"] == margin_channels
                    and row["minimum_amplitude_snr"]
                    == minimum_amplitude_snr
                ]
                if not selected:
                    continue
                passed_total = sum(
                    row["masked_fraction"]
                    <= arguments.maximum_mask_fraction
                    and row["estimated_leakage_fraction"]
                    <= arguments.maximum_leakage_fraction
                    and row["all_tooth_leakage_metrics_finite"]
                    and row["maximum_tooth_leakage_fraction"]
                    <= arguments.maximum_tooth_leakage_fraction
                    for row in selected
                )
                passed_each = sum(
                    row["masked_fraction"]
                    <= arguments.maximum_mask_fraction
                    and row["maximum_tooth_leakage_fraction"]
                    <= arguments.maximum_tooth_leakage_fraction
                    for row in selected
                )
                print(
                    f"fwhm={half_width_fwhm:g},"
                    f"margin={margin_channels:g},"
                    f"minimum_snr={minimum_amplitude_snr:g}: "
                    f"pass_total={passed_total}/{len(selected)}, "
                    f"pass_each={passed_each}/{len(selected)}, "
                    "median_teeth="
                    f"{np.median([row['masked_tooth_count'] for row in selected]):.1f}, "
                    "median_mask="
                    f"{np.median([row['masked_fraction'] for row in selected]):.6f}, "
                    "median_leakage="
                    f"{np.median([row['estimated_leakage_fraction'] for row in selected]):.6f}"
                )
    print(f"saved_csv={arguments.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
