#!/usr/bin/env python3
"""Select local-comb candidates after testing their 15/30 s tracks."""

from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml

from compare_stable_local_tracks import correlation
from discover_local_comb_guides import load_tracker_parameters
from track_dynamic_waterfall import track_comb


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--candidate-pairs", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--scan", type=int, required=True)
    parser.add_argument("--polar", choices=["XX", "YY"], required=True)
    parser.add_argument("--left-waterfall-variant", required=True)
    parser.add_argument("--right-waterfall-variant", required=True)
    parser.add_argument(
        "--maximum-disagreement-channels", type=float, default=4.0
    )
    parser.add_argument(
        "--minimum-agreement-fraction", type=float, default=0.8
    )
    parser.add_argument("--minimum-matched-guides", type=int, default=3)
    parser.add_argument("--require-both-individually-accepted", action="store_true")
    parser.add_argument("--require-spacing-regular", action="store_true")
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def load_waterfall(
    tables_dir: Path,
    scan: int,
    beam: str,
    polar: str,
    variant: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = tables_dir / (
        f"dynamic_mask_scan{scan}_{beam}_{polar}_{variant}_waterfall.npz"
    )
    with np.load(path, allow_pickle=False) as product:
        return (
            product["frequency_mhz"].astype(np.float64),
            product["binned_mjd"].astype(np.float64),
            product["smoothed_standardized"].astype(np.float64),
        )


def track_metrics(
    left_frequency_mhz: np.ndarray,
    left_mjd: np.ndarray,
    left_waterfall: np.ndarray,
    right_frequency_mhz: np.ndarray,
    right_mjd: np.ndarray,
    right_waterfall: np.ndarray,
    guides_mhz: np.ndarray,
    base_parameters,
    maximum_disagreement_channels: float,
) -> dict[str, float]:
    left_parameters = replace(
        base_parameters,
        guide_centers_mhz=guides_mhz,
        teeth_used_for_score=min(
            base_parameters.teeth_used_for_score, len(guides_mhz)
        ),
        minimum_supported_teeth=min(
            base_parameters.minimum_supported_teeth, len(guides_mhz)
        ),
    )
    left = track_comb(
        left_frequency_mhz, left_waterfall, left_parameters
    )
    right = track_comb(
        right_frequency_mhz, right_waterfall, left_parameters
    )
    overlap = (
        (right_mjd >= left_mjd[0]) & (right_mjd <= left_mjd[-1])
    )
    if not np.any(overlap):
        raise ValueError("The two Waterfalls have no overlapping time bins")
    comparison_shift = np.interp(
        right_mjd[overlap],
        left_mjd,
        left["common_shift_mhz"],
    )
    right_shift = right["common_shift_mhz"][overlap]
    delta_frequency_mhz = abs(
        float(np.median(np.diff(right_frequency_mhz)))
    )
    difference_channels = (
        np.abs(right_shift - comparison_shift) / delta_frequency_mhz
    )
    return {
        "left_reliable_fraction": float(np.mean(left["reliable"])),
        "right_reliable_fraction": float(np.mean(right["reliable"])),
        "shift_correlation": correlation(comparison_shift, right_shift),
        "difference_channels_median": float(
            np.median(difference_channels)
        ),
        "difference_channels_p95": float(
            np.percentile(difference_channels, 95)
        ),
        "difference_channels_maximum": float(
            np.max(difference_channels)
        ),
        "agreement_within_2_channels_fraction": float(
            np.mean(difference_channels <= 2)
        ),
        "agreement_within_3_channels_fraction": float(
            np.mean(difference_channels <= 3)
        ),
        "agreement_within_4_channels_fraction": float(
            np.mean(
                difference_channels <= maximum_disagreement_channels
            )
        ),
    }


def finite_or(value: str | float, fallback: float) -> float:
    parsed = float(value)
    return parsed if np.isfinite(parsed) else fallback


def ranking_key(row: dict[str, str | float | bool]) -> tuple:
    """Prefer accepted four-tooth tracks, then stronger agreement."""
    return (
        bool(row["track_consistency_accepted"]),
        int(row["matched_count"]),
        float(row["agreement_within_4_channels_fraction"]),
        float(row["agreement_within_3_channels_fraction"]),
        finite_or(row["shift_correlation"], -2.0),
        min(
            float(row["left_dynamic_path_z"]),
            float(row["right_dynamic_path_z"]),
        ),
        float(row["combined_selection_score"]),
        -float(row["maximum_difference_khz"]),
    )


def main() -> int:
    arguments = parse_arguments()
    with arguments.config.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    base_parameters = load_tracker_parameters(config["tracker"])
    with arguments.candidate_pairs.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        input_rows = list(csv.DictReader(handle))

    eligible = [
        row
        for row in input_rows
        if int(row["matched_count"]) >= arguments.minimum_matched_guides
        and (
            not arguments.require_both_individually_accepted
            or row["both_individually_accepted"].lower() == "true"
        )
        and (
            not arguments.require_spacing_regular
            or row["both_spacing_regular"].lower() == "true"
        )
    ]
    if not eligible:
        raise RuntimeError("No candidate pairs meet the pre-track requirements")

    output_rows: list[dict[str, str | float | bool]] = []
    for beam in sorted({row["beam"] for row in eligible}):
        left = load_waterfall(
            arguments.tables_dir,
            arguments.scan,
            beam,
            arguments.polar,
            arguments.left_waterfall_variant,
        )
        right = load_waterfall(
            arguments.tables_dir,
            arguments.scan,
            beam,
            arguments.polar,
            arguments.right_waterfall_variant,
        )
        beam_rows = []
        for row in (item for item in eligible if item["beam"] == beam):
            guides_mhz = np.asarray(
                [
                    float(value)
                    for value in row["matched_mean_guides_mhz"].split(";")
                ],
                dtype=np.float64,
            )
            metrics = track_metrics(
                *left,
                *right,
                guides_mhz,
                base_parameters,
                arguments.maximum_disagreement_channels,
            )
            result = {
                **row,
                **metrics,
                "track_consistency_accepted": (
                    metrics["agreement_within_4_channels_fraction"]
                    >= arguments.minimum_agreement_fraction
                ),
            }
            beam_rows.append(result)
        beam_rows.sort(key=ranking_key, reverse=True)
        for track_rank, row in enumerate(beam_rows, start=1):
            row["track_rank"] = track_rank
            row["selected"] = track_rank == 1
            output_rows.append(row)
        selected = beam_rows[0]
        print(
            f"{beam}: tested={len(beam_rows)}, "
            f"selected_pair_rank={selected['pair_rank']}, "
            f"guides={selected['matched_count']}, "
            "agreement4="
            f"{selected['agreement_within_4_channels_fraction']:.6f}, "
            f"accepted={selected['track_consistency_accepted']}"
        )

    fieldnames = [
        "track_rank",
        "selected",
        *[
            key
            for key in output_rows[0]
            if key not in {"track_rank", "selected"}
        ],
    ]
    arguments.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output_csv.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"saved_csv={arguments.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
