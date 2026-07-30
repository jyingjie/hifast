#!/usr/bin/env python3
"""Evaluate stable local-comb masks for every accepted beam."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--track-consistency-csv", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--reports-dir", type=Path, required=True)
    parser.add_argument("--scan", type=int, default=43)
    parser.add_argument("--polar", choices=["XX", "YY"], default="YY")
    parser.add_argument("--primary-label", required=True)
    parser.add_argument("--comparison-label", required=True)
    parser.add_argument("--output-label", required=True)
    parser.add_argument("--output-summary-csv", type=Path, required=True)
    parser.add_argument("--maximum-leakage-fraction", type=float, default=0.10)
    parser.add_argument("--maximum-mask-fraction", type=float, default=0.15)
    parser.add_argument(
        "--conservative-half-width-fwhm", type=float, default=0.80
    )
    parser.add_argument(
        "--conservative-margin-channels", type=float, default=0.50
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    arguments = parse_arguments()
    evaluator = Path(__file__).with_name("evaluate_local_group_mask.py")
    consistency_rows = [
        row
        for row in read_csv(arguments.track_consistency_csv)
        if row["track_consistency_accepted"].lower() == "true"
    ]
    output_rows = []
    for consistency in consistency_rows:
        beam = consistency["beam"]
        stem = (
            f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
            f"{arguments.output_label}"
        )
        profile_csv = arguments.tables_dir / f"{stem}_profile_fits.csv"
        metrics_csv = arguments.tables_dir / f"{stem}_mask_metrics.csv"
        tooth_csv = arguments.tables_dir / f"{stem}_tooth_metrics.csv"
        command = [
            sys.executable,
            str(evaluator),
            "--waterfall",
            str(
                arguments.tables_dir
                / (
                    f"dynamic_mask_scan{arguments.scan}_{beam}_"
                    f"{arguments.polar}_"
                    "pre_fc_30s_waterfall.npz"
                )
            ),
            "--primary-tracks",
            str(
                arguments.tables_dir
                / (
                    f"dynamic_mask_scan{arguments.scan}_{beam}_"
                    f"{arguments.polar}_"
                    f"{arguments.primary_label}_tracks.npz"
                )
            ),
            "--comparison-tracks",
            str(
                arguments.tables_dir
                / (
                    f"dynamic_mask_scan{arguments.scan}_{beam}_"
                    f"{arguments.polar}_"
                    f"{arguments.comparison_label}_tracks.npz"
                )
            ),
            "--sample-id",
            (
                f"scan{arguments.scan}_{beam}_{arguments.polar}_"
                f"{arguments.output_label}"
            ),
            "--group-label",
            "稳定主组",
            "--conservative-half-width-fwhm",
            str(arguments.conservative_half_width_fwhm),
            "--conservative-margin-channels",
            str(arguments.conservative_margin_channels),
            "--output-masks",
            str(arguments.tables_dir / f"{stem}_candidate_masks.npz"),
            "--profile-csv",
            str(profile_csv),
            "--metrics-csv",
            str(metrics_csv),
            "--tooth-metrics-csv",
            str(tooth_csv),
            "--report",
            str(
                arguments.reports_dir
                / (
                    f"STAGE_2_SCAN{arguments.scan}_{beam}_"
                    f"{arguments.polar}_"
                    f"{arguments.output_label}.md"
                )
            ),
            "--figure",
            str(arguments.figures_dir / f"{stem}_evaluation.png"),
        ]
        print(f"evaluate={beam}", flush=True)
        subprocess.run(command, check=True)
        profile_rows = read_csv(profile_csv)
        metric_rows = {
            row["candidate"]: row for row in read_csv(metrics_csv)
        }
        tooth_rows = [
            row
            for row in read_csv(tooth_csv)
            if row["candidate"] == "conservative_all_times"
        ]
        conservative = metric_rows["conservative_all_times"]
        snr = np.asarray(
            [float(row["amplitude_snr"]) for row in profile_rows]
        )
        widths = np.asarray(
            [float(row["fwhm_khz"]) for row in profile_rows]
        )
        leakage = float(conservative["estimated_leakage_fraction"])
        mask_fraction = float(conservative["masked_fraction"])
        tooth_leakages = np.asarray(
            [
                float(row["estimated_leakage_fraction"])
                for row in tooth_rows
            ],
            dtype=np.float64,
        )
        evaluation_complete = bool(np.all(np.isfinite(tooth_leakages)))
        maximum_tooth_leakage = float(np.nanmax(tooth_leakages))
        passes_total_targets = (
            leakage <= arguments.maximum_leakage_fraction
            and mask_fraction <= arguments.maximum_mask_fraction
        )
        output_rows.append(
            {
                "beam": beam,
                "guide_count": consistency["matched_guide_count"],
                "track_agreement_within_4_channels_fraction": consistency[
                    "agreement_within_4_channels_fraction"
                ],
                "track_shift_correlation": consistency["shift_correlation"],
                "minimum_amplitude_snr": float(np.min(snr)),
                "median_amplitude_snr": float(np.median(snr)),
                "minimum_fwhm_khz": float(np.min(widths)),
                "maximum_fwhm_khz": float(np.max(widths)),
                "conservative_mask_fraction": mask_fraction,
                "conservative_estimated_leakage_fraction": leakage,
                "maximum_tooth_leakage_fraction": maximum_tooth_leakage,
                "channel_time_loss_p95": conservative[
                    "channel_time_loss_p95"
                ],
                "channel_time_loss_max": conservative[
                    "channel_time_loss_max"
                ],
                "all_tooth_leakage_metrics_finite": evaluation_complete,
                "passes_total_leakage_and_mask_targets": (
                    passes_total_targets
                ),
                "passes_each_tooth_leakage_target": (
                    evaluation_complete
                    and maximum_tooth_leakage
                    <= arguments.maximum_leakage_fraction
                ),
                "passes_complete_mask_evaluation": (
                    evaluation_complete and passes_total_targets
                ),
            }
        )
    if not output_rows:
        raise RuntimeError("No track-consistent beams were found")
    arguments.output_summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output_summary_csv.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    for row in output_rows:
        print(
            f"{row['beam']}: mask="
            f"{float(row['conservative_mask_fraction']):.6f}, "
            "leakage="
            f"{float(row['conservative_estimated_leakage_fraction']):.6f}, "
            "maximum_tooth_leakage="
            f"{float(row['maximum_tooth_leakage_fraction']):.6f}, "
            "passes="
            f"{row['passes_total_leakage_and_mask_targets']}"
        )
    print(f"saved_summary_csv={arguments.output_summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
