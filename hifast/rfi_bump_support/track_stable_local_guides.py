#!/usr/bin/env python3
"""Track cross-time-bin stable local comb guides in multiple beams."""

from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from discover_local_comb_guides import load_tracker_parameters
from track_dynamic_waterfall import track_comb


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--candidate-pairs", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--scan", type=int, default=43)
    parser.add_argument("--polar", choices=["XX", "YY"], default="YY")
    parser.add_argument("--waterfall-variant", required=True)
    parser.add_argument("--output-label", required=True)
    parser.add_argument("--minimum-matched-guides", type=int, default=3)
    parser.add_argument(
        "--require-both-individually-accepted", action="store_true"
    )
    parser.add_argument("--require-spacing-regular", action="store_true")
    return parser.parse_args()


def selected_pairs(
    path: Path,
    minimum_matched_guides: int,
    require_both_accepted: bool,
    require_spacing_regular: bool,
) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle)
            if row["selected"].lower() == "true"
            and int(row["matched_count"]) >= minimum_matched_guides
            and (
                not require_both_accepted
                or row["both_individually_accepted"].lower() == "true"
            )
            and (
                not require_spacing_regular
                or row["both_spacing_regular"].lower() == "true"
            )
        ]
    return rows


def plot_track(
    path: Path,
    scan: int,
    beam: str,
    polar: str,
    frequency_mhz: np.ndarray,
    elapsed_minutes: np.ndarray,
    waterfall: np.ndarray,
    guides_mhz: np.ndarray,
    result: dict[str, np.ndarray],
) -> None:
    figure, axes = plt.subplots(
        2,
        1,
        figsize=(14, 8),
        gridspec_kw={"height_ratios": [3.0, 1.0]},
        constrained_layout=True,
    )
    limit = max(2.0, float(np.nanpercentile(np.abs(waterfall), 99.0)))
    image = axes[0].imshow(
        waterfall,
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
        vmin=-limit,
        vmax=limit,
    )
    for tooth_index in range(len(guides_mhz)):
        axes[0].plot(
            result["common_track_mhz"][:, tooth_index],
            elapsed_minutes,
            color="black",
            linewidth=1.0,
        )
    axes[0].set(
        ylabel="time from start (min)",
        title=f"scan {scan} {beam}/{polar}: stable local guides",
    )
    figure.colorbar(
        image,
        ax=axes[0],
        pad=0.01,
        label="smoothed response / robust RMS",
    )
    axes[1].plot(
        elapsed_minutes,
        result["common_shift_mhz"] * 1000.0,
        color="black",
    )
    axes[1].scatter(
        elapsed_minutes[result["reliable"]],
        result["common_shift_mhz"][result["reliable"]] * 1000.0,
        s=10,
        color="tab:green",
        label="reliable time bin",
    )
    axes[1].set(
        xlabel="time from start (min)",
        ylabel="common shift (kHz)",
    )
    axes[1].grid(alpha=0.2)
    axes[1].legend()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> int:
    arguments = parse_arguments()
    with arguments.config.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    base_parameters = load_tracker_parameters(config["tracker"])
    rows = selected_pairs(
        arguments.candidate_pairs,
        arguments.minimum_matched_guides,
        arguments.require_both_individually_accepted,
        arguments.require_spacing_regular,
    )
    if not rows:
        raise RuntimeError("No selected candidate pairs meet the requirements")

    for row in rows:
        beam = row["beam"]
        guides_mhz = np.asarray(
            [
                float(value)
                for value in row["matched_mean_guides_mhz"].split(";")
            ],
            dtype=np.float64,
        )
        waterfall_path = arguments.tables_dir / (
            f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
            f"{arguments.waterfall_variant}_waterfall.npz"
        )
        with np.load(waterfall_path, allow_pickle=False) as product:
            frequency_mhz = product["frequency_mhz"].astype(np.float64)
            waterfall = product[
                "smoothed_standardized"
            ].astype(np.float64)
            binned_mjd = product["binned_mjd"].astype(np.float64)
            elapsed_minutes = product[
                "binned_elapsed_minutes"
            ].astype(np.float64)
        parameters = replace(
            base_parameters,
            guide_centers_mhz=guides_mhz,
            teeth_used_for_score=min(
                base_parameters.teeth_used_for_score, len(guides_mhz)
            ),
            minimum_supported_teeth=min(
                base_parameters.minimum_supported_teeth, len(guides_mhz)
            ),
        )
        result = track_comb(frequency_mhz, waterfall, parameters)
        stem = (
            f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
            f"{arguments.output_label}"
        )
        output_npz = arguments.tables_dir / f"{stem}_tracks.npz"
        np.savez_compressed(
            output_npz,
            beam=beam,
            accepted_pair=(
                row["both_individually_accepted"].lower() == "true"
            ),
            matched_guide_count=int(row["matched_count"]),
            maximum_guide_difference_khz=float(
                row["maximum_difference_khz"]
            ),
            frequency_mhz=frequency_mhz,
            binned_mjd=binned_mjd,
            binned_elapsed_minutes=elapsed_minutes,
            guide_centers_mhz=guides_mhz,
            independent_common_shift_mhz=result["common_shift_mhz"],
            independent_common_track_mhz=result["common_track_mhz"],
            independent_reliable=result["reliable"],
            independent_path_confidence=result["path_confidence"],
            independent_support_count=result["support_count"],
        )
        plot_track(
            arguments.figures_dir / f"{stem}.png",
            arguments.scan,
            beam,
            arguments.polar,
            frequency_mhz,
            elapsed_minutes,
            waterfall,
            guides_mhz,
            result,
        )
        print(
            f"{beam}: guides={len(guides_mhz)}, "
            f"reliable={np.mean(result['reliable']):.6f}, "
            f"saved={output_npz}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
