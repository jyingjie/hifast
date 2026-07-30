#!/usr/bin/env python3
"""Track a commonly drifting comb in a floating-point Waterfall array."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy.ndimage import median_filter


@dataclass(frozen=True)
class TrackerParameters:
    guide_centers_mhz: np.ndarray
    maximum_common_shift_mhz: float
    maximum_jump_channels: int
    transition_penalty: float
    teeth_used_for_score: int
    response_clip_min: float
    response_clip_max: float
    local_refine_half_width_channels: int
    common_refine_time_median_rows: int
    support_response_threshold: float
    minimum_supported_teeth: int
    minimum_path_confidence: float


def robust_scale(values: np.ndarray, axis: int | None = None) -> np.ndarray:
    """Return a MAD-based scale with a finite nonzero fallback."""
    center = np.nanmedian(values, axis=axis, keepdims=True)
    scale = 1.4826 * np.nanmedian(
        np.abs(values - center), axis=axis, keepdims=True
    )
    standard_deviation = np.nanstd(values, axis=axis, keepdims=True)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, standard_deviation)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, 1.0)
    if axis is None:
        return np.asarray(scale).reshape(())
    return np.squeeze(scale, axis=axis)


def nearest_channels(
    frequency_mhz: np.ndarray, centers_mhz: np.ndarray
) -> np.ndarray:
    """Return the nearest channel for each requested frequency."""
    if frequency_mhz.ndim != 1 or len(frequency_mhz) < 2:
        raise ValueError("frequency_mhz must be a one-dimensional axis")
    if not np.all(np.diff(frequency_mhz) > 0):
        raise ValueError("frequency_mhz must be strictly increasing")
    right = np.searchsorted(frequency_mhz, centers_mhz)
    right = np.clip(right, 1, len(frequency_mhz) - 1)
    left = right - 1
    choose_right = (
        np.abs(frequency_mhz[right] - centers_mhz)
        < np.abs(frequency_mhz[left] - centers_mhz)
    )
    return np.where(choose_right, right, left).astype(np.int64)


def build_response_cube(
    waterfall: np.ndarray,
    guide_channels: np.ndarray,
    state_offsets: np.ndarray,
) -> np.ndarray:
    """Sample every comb tooth at every possible common channel shift."""
    channels = guide_channels[None, :] + state_offsets[:, None]
    if np.any(channels < 0) or np.any(channels >= waterfall.shape[1]):
        raise ValueError("The common-shift search reaches a frequency edge")
    return waterfall[:, channels]


def emission_scores(
    response_cube: np.ndarray,
    teeth_used_for_score: int,
    clip_min: float,
    clip_max: float,
) -> np.ndarray:
    """Calculate a robust common-comb score for each time and shift."""
    tooth_count = response_cube.shape[2]
    if not 1 <= teeth_used_for_score <= tooth_count:
        raise ValueError("teeth_used_for_score is outside the tooth count")
    clipped = np.clip(response_cube, clip_min, clip_max)
    sorted_responses = np.sort(clipped, axis=2)
    selected = sorted_responses[:, :, -teeth_used_for_score:]
    raw_score = np.nanmean(selected, axis=2)
    row_center = np.nanmedian(raw_score, axis=1, keepdims=True)
    row_scale = robust_scale(raw_score, axis=1)[:, None]
    return (raw_score - row_center) / row_scale


def viterbi_path(
    emission: np.ndarray,
    state_offsets: np.ndarray,
    maximum_jump_channels: int,
    transition_penalty: float,
) -> np.ndarray:
    """Find the maximum-score continuous common-shift path."""
    time_count, state_count = emission.shape
    cumulative = np.full_like(emission, -np.inf, dtype=np.float64)
    previous = np.full((time_count, state_count), -1, dtype=np.int32)
    cumulative[0] = emission[0]

    for time_index in range(1, time_count):
        for state_index in range(state_count):
            distance = state_offsets[state_index] - state_offsets
            allowed = np.abs(distance) <= maximum_jump_channels
            candidate = cumulative[time_index - 1, allowed]
            candidate = candidate - transition_penalty * distance[allowed] ** 2
            local_previous = int(np.argmax(candidate))
            previous_states = np.flatnonzero(allowed)
            best_previous = int(previous_states[local_previous])
            cumulative[time_index, state_index] = (
                emission[time_index, state_index]
                + candidate[local_previous]
            )
            previous[time_index, state_index] = best_previous

    path_state = np.empty(time_count, dtype=np.int64)
    path_state[-1] = int(np.argmax(cumulative[-1]))
    for time_index in range(time_count - 1, 0, -1):
        path_state[time_index - 1] = previous[
            time_index, path_state[time_index]
        ]
    return path_state


def refine_common_shift(
    waterfall: np.ndarray,
    guide_channels: np.ndarray,
    initial_shift_channels: np.ndarray,
    half_width_channels: int,
    time_median_rows: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Refine tooth peaks locally and return their median common shift."""
    time_count = waterfall.shape[0]
    tooth_count = len(guide_channels)
    local_offsets = np.arange(-half_width_channels, half_width_channels + 1)
    local_peak_offsets = np.zeros((time_count, tooth_count), dtype=np.int64)
    local_peak_responses = np.full(
        (time_count, tooth_count), np.nan, dtype=np.float64
    )

    for time_index in range(time_count):
        expected = guide_channels + initial_shift_channels[time_index]
        for tooth_index, expected_channel in enumerate(expected):
            channels = expected_channel + local_offsets
            valid = (channels >= 0) & (channels < waterfall.shape[1])
            values = waterfall[time_index, channels[valid]]
            if not np.any(np.isfinite(values)):
                continue
            local_index = int(np.nanargmax(values))
            local_peak_offsets[time_index, tooth_index] = local_offsets[
                valid
            ][local_index]
            local_peak_responses[time_index, tooth_index] = values[local_index]

    median_local_offset = np.rint(
        np.nanmedian(local_peak_offsets, axis=1)
    ).astype(np.int64)
    refined = initial_shift_channels + median_local_offset
    if time_median_rows > 1:
        if time_median_rows % 2 == 0:
            raise ValueError("common_refine_time_median_rows must be odd")
        refined = median_filter(
            refined, size=time_median_rows, mode="nearest"
        ).astype(np.int64)
    return refined, local_peak_offsets, local_peak_responses


def track_comb(
    frequency_mhz: np.ndarray,
    waterfall: np.ndarray,
    parameters: TrackerParameters,
) -> dict[str, np.ndarray]:
    """Track the comb and return common and per-tooth diagnostics."""
    delta_frequency_mhz = float(np.median(np.diff(frequency_mhz)))
    guide_channels = nearest_channels(
        frequency_mhz, parameters.guide_centers_mhz
    )
    maximum_shift_channels = int(
        np.floor(
            parameters.maximum_common_shift_mhz
            / abs(delta_frequency_mhz)
        )
    )
    edge_limit = min(
        int(np.min(guide_channels)),
        int(len(frequency_mhz) - 1 - np.max(guide_channels)),
    )
    maximum_shift_channels = min(maximum_shift_channels, edge_limit)
    if maximum_shift_channels < 1:
        raise ValueError("No common-shift search range is available")
    state_offsets = np.arange(
        -maximum_shift_channels,
        maximum_shift_channels + 1,
        dtype=np.int64,
    )

    response_cube = build_response_cube(
        waterfall, guide_channels, state_offsets
    )
    emission = emission_scores(
        response_cube,
        parameters.teeth_used_for_score,
        parameters.response_clip_min,
        parameters.response_clip_max,
    )
    path_state = viterbi_path(
        emission,
        state_offsets,
        parameters.maximum_jump_channels,
        parameters.transition_penalty,
    )
    initial_shift_channels = state_offsets[path_state]
    (
        refined_shift_channels,
        local_peak_offsets,
        local_peak_responses,
    ) = refine_common_shift(
        waterfall,
        guide_channels,
        initial_shift_channels,
        parameters.local_refine_half_width_channels,
        parameters.common_refine_time_median_rows,
    )

    time_indices = np.arange(waterfall.shape[0])
    selected_responses = response_cube[time_indices, path_state]
    selected_emission = emission[time_indices, path_state]
    row_center = np.nanmedian(emission, axis=1)
    row_scale = robust_scale(emission, axis=1)
    path_confidence = (selected_emission - row_center) / row_scale
    support_count = np.count_nonzero(
        local_peak_responses >= parameters.support_response_threshold,
        axis=1,
    )
    reliable = (
        (support_count >= parameters.minimum_supported_teeth)
        & (path_confidence >= parameters.minimum_path_confidence)
    )

    common_shift_mhz = refined_shift_channels * delta_frequency_mhz
    common_track_mhz = (
        parameters.guide_centers_mhz[None, :]
        + common_shift_mhz[:, None]
    )
    local_track_channels = (
        guide_channels[None, :]
        + initial_shift_channels[:, None]
        + local_peak_offsets
    )
    local_track_channels = np.clip(
        local_track_channels, 0, len(frequency_mhz) - 1
    )
    local_track_mhz = frequency_mhz[local_track_channels]

    return {
        "guide_channels": guide_channels,
        "state_offsets": state_offsets,
        "emission": emission,
        "initial_shift_channels": initial_shift_channels,
        "common_shift_channels": refined_shift_channels,
        "common_shift_mhz": common_shift_mhz,
        "common_track_mhz": common_track_mhz,
        "local_track_mhz": local_track_mhz,
        "selected_tooth_responses": selected_responses,
        "local_peak_responses": local_peak_responses,
        "support_count": support_count,
        "path_confidence": path_confidence,
        "reliable": reliable,
        "delta_frequency_mhz": np.asarray(delta_frequency_mhz),
    }


def load_parameters(
    path: Path, group: str | None = None
) -> tuple[TrackerParameters, dict]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    values = config["comb_tracker"]
    guide_centers = values["guide_centers_mhz"]
    if group is not None:
        groups = values.get("groups", {})
        if group not in groups:
            available = ", ".join(sorted(groups)) or "none"
            raise KeyError(
                f"Unknown comb group {group!r}; available groups: {available}"
            )
        guide_centers = groups[group]
    parameters = TrackerParameters(
        guide_centers_mhz=np.asarray(
            guide_centers, dtype=np.float64
        ),
        maximum_common_shift_mhz=float(
            values["maximum_common_shift_mhz"]
        ),
        maximum_jump_channels=int(
            values["maximum_jump_channels_per_time_bin"]
        ),
        transition_penalty=float(
            values["transition_penalty_per_channel_squared"]
        ),
        teeth_used_for_score=int(values["teeth_used_for_score"]),
        response_clip_min=float(values["response_clip"][0]),
        response_clip_max=float(values["response_clip"][1]),
        local_refine_half_width_channels=int(
            values["local_refine_half_width_channels"]
        ),
        common_refine_time_median_rows=int(
            values["common_refine_time_median_rows"]
        ),
        support_response_threshold=float(
            values["support_response_threshold"]
        ),
        minimum_supported_teeth=int(values["minimum_supported_teeth"]),
        minimum_path_confidence=float(
            values["minimum_path_confidence"]
        ),
    )
    return parameters, config


def write_track_csv(
    path: Path,
    elapsed_minutes: np.ndarray,
    result: dict[str, np.ndarray],
) -> None:
    tooth_count = result["common_track_mhz"].shape[1]
    fieldnames = [
        "time_bin",
        "elapsed_minutes",
        "common_shift_channels",
        "common_shift_khz",
        "path_confidence",
        "support_count",
        "reliable",
    ]
    for tooth_index in range(tooth_count):
        number = tooth_index + 1
        fieldnames.extend(
            [
                f"tooth_{number}_common_frequency_mhz",
                f"tooth_{number}_local_frequency_mhz",
                f"tooth_{number}_response",
            ]
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for time_index, elapsed in enumerate(elapsed_minutes):
            row: dict[str, object] = {
                "time_bin": time_index,
                "elapsed_minutes": f"{elapsed:.9f}",
                "common_shift_channels": int(
                    result["common_shift_channels"][time_index]
                ),
                "common_shift_khz": (
                    f"{result['common_shift_mhz'][time_index] * 1000:.6f}"
                ),
                "path_confidence": (
                    f"{result['path_confidence'][time_index]:.6f}"
                ),
                "support_count": int(result["support_count"][time_index]),
                "reliable": bool(result["reliable"][time_index]),
            }
            for tooth_index in range(tooth_count):
                number = tooth_index + 1
                row[f"tooth_{number}_common_frequency_mhz"] = (
                    f"{result['common_track_mhz'][time_index, tooth_index]:.9f}"
                )
                row[f"tooth_{number}_local_frequency_mhz"] = (
                    f"{result['local_track_mhz'][time_index, tooth_index]:.9f}"
                )
                row[f"tooth_{number}_response"] = (
                    f"{result['local_peak_responses'][time_index, tooth_index]:.6f}"
                )
            writer.writerow(row)


def write_review_template(
    path: Path,
    sample_id: str,
    elapsed_minutes: np.ndarray,
    result: dict[str, np.ndarray],
    stride: int,
) -> None:
    fieldnames = [
        "sample_id",
        "time_bin",
        "elapsed_minutes",
        "tooth_number",
        "automatic_frequency_mhz",
        "reference_frequency_mhz",
        "reference_left_mhz",
        "reference_right_mhz",
        "review_status",
        "notes",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for time_index in range(0, len(elapsed_minutes), stride):
            for tooth_index, frequency in enumerate(
                result["common_track_mhz"][time_index], start=1
            ):
                writer.writerow(
                    {
                        "sample_id": sample_id,
                        "time_bin": time_index,
                        "elapsed_minutes": f"{elapsed_minutes[time_index]:.9f}",
                        "tooth_number": tooth_index,
                        "automatic_frequency_mhz": f"{frequency:.9f}",
                        "reference_frequency_mhz": "",
                        "reference_left_mhz": "",
                        "reference_right_mhz": "",
                        "review_status": "pending",
                        "notes": "",
                    }
                )


def plot_result(
    path: Path,
    frequency_mhz: np.ndarray,
    elapsed_minutes: np.ndarray,
    waterfall: np.ndarray,
    result: dict[str, np.ndarray],
    sample_id: str,
) -> None:
    figure, axes = plt.subplots(
        3,
        1,
        figsize=(15, 11),
        sharex=False,
        gridspec_kw={"height_ratios": [2.1, 0.8, 0.9]},
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
    for tooth_index in range(result["common_track_mhz"].shape[1]):
        axes[0].plot(
            result["common_track_mhz"][:, tooth_index],
            elapsed_minutes,
            color="black",
            lw=1.1,
        )
        axes[0].scatter(
            result["local_track_mhz"][:, tooth_index],
            elapsed_minutes,
            c=np.where(result["reliable"], "lime", "gold"),
            s=5,
            alpha=0.65,
            edgecolors="none",
        )
    axes[0].set_xlabel("topocentric frequency (MHz)")
    axes[0].set_ylabel("time from start (min)")
    axes[0].set_title(
        "Floating-point Waterfall with common-drift tracks; "
        "black=model, green=reliable local peak, gold=low confidence"
    )
    figure.colorbar(
        image,
        ax=axes[0],
        pad=0.01,
        label="smoothed response / robust RMS",
    )

    axes[1].plot(
        elapsed_minutes,
        result["common_shift_mhz"] * 1000,
        color="tab:blue",
        lw=1.2,
    )
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].set_ylabel("common shift (kHz)")
    axes[1].grid(alpha=0.2)

    axes[2].plot(
        elapsed_minutes,
        result["path_confidence"],
        color="tab:purple",
        lw=1.0,
        label="path confidence",
    )
    axes[2].step(
        elapsed_minutes,
        result["support_count"],
        color="tab:green",
        lw=1.0,
        where="mid",
        label="supported teeth",
    )
    axes[2].scatter(
        elapsed_minutes[~result["reliable"]],
        result["path_confidence"][~result["reliable"]],
        color="tab:orange",
        s=12,
        label="low-confidence row",
    )
    axes[2].axhline(1.0, color="0.3", lw=0.6, linestyle="--")
    axes[2].set_xlabel("time from start (min)")
    axes[2].set_ylabel("confidence / count")
    axes[2].grid(alpha=0.2)
    axes[2].legend(fontsize=8, ncol=3)

    figure.suptitle(sample_id)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=170)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--group",
        help="Optional named comb group from comb_tracker.groups",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-npz", type=Path, required=True)
    parser.add_argument("--tracks-csv", type=Path, required=True)
    parser.add_argument("--review-template-csv", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()

    parameters, config = load_parameters(args.config, args.group)
    with np.load(args.input, allow_pickle=False) as product:
        frequency_mhz = product["frequency_mhz"].astype(np.float64)
        elapsed_minutes = product["binned_elapsed_minutes"].astype(np.float64)
        waterfall = product["smoothed_standardized"].astype(np.float64)
        source_polarization = str(product["polarization"])
        samples_per_bin = int(product["samples_per_bin"])
        sample_seconds = float(product["sample_seconds"])

    result = track_comb(frequency_mhz, waterfall, parameters)
    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        sample_id=args.sample_id,
        comb_group=args.group or "default",
        source_waterfall_npz=str(args.input),
        source_polarization=source_polarization,
        frequency_mhz=frequency_mhz,
        binned_elapsed_minutes=elapsed_minutes,
        samples_per_bin=samples_per_bin,
        sample_seconds=sample_seconds,
        guide_centers_mhz=parameters.guide_centers_mhz,
        **result,
    )
    write_track_csv(args.tracks_csv, elapsed_minutes, result)
    write_review_template(
        args.review_template_csv,
        args.sample_id,
        elapsed_minutes,
        result,
        int(config["review"]["proposal_stride_time_bins"]),
    )
    plot_result(
        args.figure,
        frequency_mhz,
        elapsed_minutes,
        waterfall,
        result,
        args.sample_id,
    )

    shift_khz = result["common_shift_mhz"] * 1000
    print(f"sample_id={args.sample_id}")
    print(f"time_bins={len(elapsed_minutes)}")
    print(f"time_bin_seconds={samples_per_bin * sample_seconds:.6f}")
    print(
        "common_shift_khz_min_median_max="
        f"{np.min(shift_khz):.6f},"
        f"{np.median(shift_khz):.6f},"
        f"{np.max(shift_khz):.6f}"
    )
    print(
        "common_shift_peak_to_peak_khz="
        f"{np.ptp(shift_khz):.6f}"
    )
    print(
        "reliable_time_fraction="
        f"{np.mean(result['reliable']):.6f}"
    )
    print(f"saved_npz={args.output_npz}")
    print(f"saved_tracks_csv={args.tracks_csv}")
    print(f"saved_review_template={args.review_template_csv}")
    print(f"saved_figure={args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
