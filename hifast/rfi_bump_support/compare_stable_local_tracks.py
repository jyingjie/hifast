#!/usr/bin/env python3
"""Compare selected local-comb tracks from two time-bin lengths."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-pairs", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--scan", type=int, default=43)
    parser.add_argument("--polar", choices=["XX", "YY"], default="YY")
    parser.add_argument("--left-label", required=True)
    parser.add_argument("--right-label", required=True)
    parser.add_argument("--maximum-disagreement-channels", type=float, default=4)
    parser.add_argument("--minimum-agreement-fraction", type=float, default=0.8)
    parser.add_argument("--require-spacing-regular", action="store_true")
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    if (
        len(left) < 3
        or np.nanstd(left) == 0
        or np.nanstd(right) == 0
    ):
        return np.nan
    return float(np.corrcoef(left, right)[0, 1])


def main() -> int:
    arguments = parse_arguments()
    with arguments.candidate_pairs.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        pairs = [
            row
            for row in csv.DictReader(handle)
            if row["selected"].lower() == "true"
            and row["both_individually_accepted"].lower() == "true"
            and int(row["matched_count"]) >= 3
            and (
                not arguments.require_spacing_regular
                or row["both_spacing_regular"].lower() == "true"
            )
        ]
    rows = []
    for pair in pairs:
        beam = pair["beam"]
        paths = [
            arguments.tables_dir
            / (
                f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
                f"{label}_tracks.npz"
            )
            for label in (arguments.left_label, arguments.right_label)
        ]
        with np.load(paths[0], allow_pickle=False) as product:
            left_mjd = product["binned_mjd"].astype(np.float64)
            left_track = product[
                "independent_common_track_mhz"
            ].astype(np.float64)
            left_shift = product[
                "independent_common_shift_mhz"
            ].astype(np.float64)
            left_reliable = product[
                "independent_reliable"
            ].astype(bool)
        with np.load(paths[1], allow_pickle=False) as product:
            right_frequency = product["frequency_mhz"].astype(np.float64)
            right_mjd = product["binned_mjd"].astype(np.float64)
            right_track = product[
                "independent_common_track_mhz"
            ].astype(np.float64)
            right_shift = product[
                "independent_common_shift_mhz"
            ].astype(np.float64)
            right_reliable = product[
                "independent_reliable"
            ].astype(bool)
        overlap = (
            (right_mjd >= left_mjd[0]) & (right_mjd <= left_mjd[-1])
        )
        comparison_track = np.column_stack(
            [
                np.interp(
                    right_mjd[overlap],
                    left_mjd,
                    left_track[:, tooth_index],
                )
                for tooth_index in range(left_track.shape[1])
            ]
        )
        difference_channels = (
            np.abs(right_track[overlap] - comparison_track)
            / abs(float(np.median(np.diff(right_frequency))))
        )
        maximum_per_time = np.max(difference_channels, axis=1)
        comparison_shift = np.interp(
            right_mjd[overlap], left_mjd, left_shift
        )
        threshold = arguments.maximum_disagreement_channels
        agreement = maximum_per_time <= threshold
        agreement_fraction = float(np.mean(agreement))
        rows.append(
            {
                "beam": beam,
                "matched_guide_count": int(pair["matched_count"]),
                "maximum_guide_difference_khz": float(
                    pair["maximum_difference_khz"]
                ),
                "left_reliable_fraction": float(np.mean(left_reliable)),
                "right_reliable_fraction": float(np.mean(right_reliable)),
                "shift_correlation": correlation(
                    comparison_shift, right_shift[overlap]
                ),
                "difference_channels_median": float(
                    np.median(maximum_per_time)
                ),
                "difference_channels_p95": float(
                    np.percentile(maximum_per_time, 95)
                ),
                "difference_channels_maximum": float(
                    np.max(maximum_per_time)
                ),
                "agreement_within_2_channels_fraction": float(
                    np.mean(maximum_per_time <= 2)
                ),
                "agreement_within_3_channels_fraction": float(
                    np.mean(maximum_per_time <= 3)
                ),
                "agreement_within_4_channels_fraction": agreement_fraction,
                "right_reliable_and_agreement_fraction": float(
                    np.mean(right_reliable[overlap] & agreement)
                ),
                "track_consistency_accepted": (
                    agreement_fraction
                    >= arguments.minimum_agreement_fraction
                ),
            }
        )
    if not rows:
        raise RuntimeError("No accepted candidate pairs were found")
    arguments.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output_csv.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(
            f"{row['beam']}: guides={row['matched_guide_count']}, "
            "agreement4="
            f"{row['agreement_within_4_channels_fraction']:.6f}, "
            f"r={row['shift_correlation']:.6f}, "
            f"accepted={row['track_consistency_accepted']}"
        )
    print(f"saved_csv={arguments.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
