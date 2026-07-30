#!/usr/bin/env python3
"""Select a comb candidate that repeats across two time-bin lengths."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left-candidates", type=Path, required=True)
    parser.add_argument("--right-candidates", type=Path, required=True)
    parser.add_argument(
        "--maximum-match-khz", type=float, default=50.0
    )
    parser.add_argument(
        "--maximum-spacing-std-mhz", type=float, default=0.075
    )
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def read_candidates(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                {
                    **raw,
                    "rank": int(raw["rank"]),
                    "tooth_count": int(raw["tooth_count"]),
                    "guides": np.asarray(
                        [
                            float(value)
                            for value in raw[
                                "guide_centers_mhz"
                            ].split(";")
                        ]
                    ),
                    "dynamic_path_z": float(raw["dynamic_path_z"]),
                    "search_edge_fraction": float(
                        raw["search_edge_fraction"]
                    ),
                    "selection_score": float(raw["selection_score"]),
                    "accepted": raw["accepted"].lower() == "true",
                    "spacing_std_mhz": float(
                        np.std(
                            np.diff(
                                [
                                    float(value)
                                    for value in raw[
                                        "guide_centers_mhz"
                                    ].split(";")
                                ]
                            )
                        )
                    ),
                }
            )
    return rows


def main() -> int:
    arguments = parse_arguments()
    left_rows = read_candidates(arguments.left_candidates)
    right_rows = read_candidates(arguments.right_candidates)
    maximum_match_mhz = arguments.maximum_match_khz / 1000.0
    pairs = []
    for left in left_rows:
        for right in right_rows:
            if left["beam"] != right["beam"]:
                continue
            distance = np.abs(
                left["guides"][:, None] - right["guides"][None, :]
            )
            left_indices, right_indices = linear_sum_assignment(distance)
            keep = distance[left_indices, right_indices] <= maximum_match_mhz
            matched_left = left_indices[keep]
            matched_right = right_indices[keep]
            differences = distance[matched_left, matched_right]
            matched_count = len(differences)
            if not matched_count:
                continue
            pairs.append(
                {
                    "beam": left["beam"],
                    "left_rank": left["rank"],
                    "right_rank": right["rank"],
                    "left_tooth_count": left["tooth_count"],
                    "right_tooth_count": right["tooth_count"],
                    "matched_count": matched_count,
                    "matched_fraction_of_smaller_set": (
                        matched_count
                        / min(
                            left["tooth_count"], right["tooth_count"]
                        )
                    ),
                    "maximum_difference_khz": (
                        float(np.max(differences) * 1000.0)
                    ),
                    "mean_difference_khz": (
                        float(np.mean(differences) * 1000.0)
                    ),
                    "left_dynamic_path_z": left["dynamic_path_z"],
                    "right_dynamic_path_z": right["dynamic_path_z"],
                    "left_search_edge_fraction": (
                        left["search_edge_fraction"]
                    ),
                    "right_search_edge_fraction": (
                        right["search_edge_fraction"]
                    ),
                    "both_individually_accepted": (
                        left["accepted"] and right["accepted"]
                    ),
                    "left_spacing_std_mhz": left["spacing_std_mhz"],
                    "right_spacing_std_mhz": right["spacing_std_mhz"],
                    "both_spacing_regular": (
                        left["spacing_std_mhz"]
                        <= arguments.maximum_spacing_std_mhz
                        and right["spacing_std_mhz"]
                        <= arguments.maximum_spacing_std_mhz
                    ),
                    "left_guides_mhz": left["guide_centers_mhz"],
                    "right_guides_mhz": right["guide_centers_mhz"],
                    "matched_mean_guides_mhz": ";".join(
                        f"{value:.9f}"
                        for value in 0.5
                        * (
                            left["guides"][matched_left]
                            + right["guides"][matched_right]
                        )
                    ),
                    "combined_selection_score": (
                        left["selection_score"]
                        + right["selection_score"]
                    ),
                }
            )
    if not pairs:
        raise RuntimeError("No candidate pairs have matching guides")

    def ranking_key(row: dict[str, object]) -> tuple:
        return (
            bool(row["both_individually_accepted"]),
            bool(row["both_spacing_regular"]),
            int(row["matched_count"]),
            float(row["matched_fraction_of_smaller_set"]),
            -float(row["maximum_difference_khz"]),
            min(
                float(row["left_dynamic_path_z"]),
                float(row["right_dynamic_path_z"]),
            ),
            float(row["combined_selection_score"]),
        )

    ranked_pairs = []
    for beam in sorted({str(row["beam"]) for row in pairs}):
        beam_pairs = [
            row for row in pairs if str(row["beam"]) == beam
        ]
        beam_pairs.sort(key=ranking_key, reverse=True)
        for pair_rank, row in enumerate(beam_pairs, start=1):
            row["pair_rank"] = pair_rank
            row["selected"] = pair_rank == 1
            ranked_pairs.append(row)
    pairs = ranked_pairs

    fieldnames = [
        "pair_rank",
        "selected",
        *[
            key
            for key in pairs[0]
            if key not in {"pair_rank", "selected"}
        ],
    ]
    arguments.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output_csv.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(pairs)

    for selected in (row for row in pairs if row["selected"]):
        print(
            f"{selected['beam']}: "
            f"selected_left_rank={selected['left_rank']}, "
            f"selected_right_rank={selected['right_rank']}, "
            f"matched_count={selected['matched_count']}, "
            "spacing_regular="
            f"{selected['both_spacing_regular']}, "
            "matched_mean_guides_mhz="
            f"{selected['matched_mean_guides_mhz']}"
        )
    print(f"saved_csv={arguments.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
