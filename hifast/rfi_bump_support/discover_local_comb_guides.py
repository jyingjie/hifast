#!/usr/bin/env python3
"""Discover beam-local comb guides and compare their drift with a reference."""

from __future__ import annotations

import argparse
import csv
import itertools
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy.signal import find_peaks

from track_dynamic_waterfall import TrackerParameters, track_comb


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--scan", type=int, default=43)
    parser.add_argument("--polar", choices=["XX", "YY"], default="YY")
    parser.add_argument("--reference-beam", default="M05")
    parser.add_argument(
        "--reference-polar", choices=["XX", "YY"], default=None
    )
    parser.add_argument("--independent-only", action="store_true")
    parser.add_argument("--reference-variant", default="pre_fc")
    parser.add_argument("--waterfall-variant", default="pre_fc")
    parser.add_argument("--output-label", default="v2_local_guides")
    parser.add_argument("--output-candidate-rank", type=int, default=1)
    parser.add_argument("--exclude-mask-product", type=Path)
    parser.add_argument(
        "--exclude-half-width-mhz", type=float, default=0.15
    )
    parser.add_argument(
        "--beams", nargs="+", default=["M02", "M07", "M13"]
    )
    parser.add_argument("--candidate-csv", type=Path, required=True)
    parser.add_argument("--summary-csv", type=Path, required=True)
    return parser.parse_args()


def robust_scale(values: np.ndarray, axis: int | None = None) -> np.ndarray:
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


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    valid = np.isfinite(left) & np.isfinite(right)
    if (
        np.count_nonzero(valid) < 3
        or np.nanstd(left[valid]) == 0
        or np.nanstd(right[valid]) == 0
    ):
        return np.nan
    return float(np.corrcoef(left[valid], right[valid])[0, 1])


def nearest_channels(
    frequency_mhz: np.ndarray, frequencies_mhz: np.ndarray
) -> np.ndarray:
    right = np.searchsorted(frequency_mhz, frequencies_mhz)
    right = np.clip(right, 1, len(frequency_mhz) - 1)
    left = right - 1
    choose_right = (
        np.abs(frequency_mhz[right] - frequencies_mhz)
        < np.abs(frequency_mhz[left] - frequencies_mhz)
    )
    return np.where(choose_right, right, left).astype(np.int64)


def load_tracker_parameters(values: dict) -> TrackerParameters:
    return TrackerParameters(
        guide_centers_mhz=np.asarray([], dtype=np.float64),
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


def build_profile(waterfall: np.ndarray) -> np.ndarray:
    profile = 0.5 * (
        np.nanmean(waterfall, axis=0)
        + np.nanmedian(waterfall, axis=0)
    )
    return profile - np.nanmedian(profile)


def profile_peaks(
    frequency_mhz: np.ndarray, profile: np.ndarray, config: dict
) -> np.ndarray:
    delta_frequency_mhz = float(np.median(np.diff(frequency_mhz)))
    minimum_distance = max(
        1,
        int(
            round(
                float(config["minimum_peak_separation_mhz"])
                / delta_frequency_mhz
            )
        ),
    )
    peaks, _ = find_peaks(
        profile,
        distance=minimum_distance,
        prominence=float(config["minimum_peak_prominence"]),
        height=float(config["minimum_peak_height"]),
    )
    edge = float(config["frequency_edge_exclusion_mhz"])
    peaks = peaks[
        (frequency_mhz[peaks] >= frequency_mhz[0] + edge)
        & (frequency_mhz[peaks] <= frequency_mhz[-1] - edge)
    ]
    maximum = int(config["maximum_profile_peaks"])
    if len(peaks) > maximum:
        peaks = peaks[np.argsort(profile[peaks])[-maximum:]]
    return np.sort(peaks)


def candidate_chains(
    frequency_mhz: np.ndarray, peaks: np.ndarray, config: dict
) -> list[tuple[np.ndarray, np.ndarray]]:
    minimum_spacing, maximum_spacing = (
        float(value) for value in config["spacing_range_mhz"]
    )
    chains = []
    for length in config["chain_lengths"]:
        for indices in itertools.combinations(peaks, int(length)):
            indices_array = np.asarray(indices, dtype=np.int64)
            guides = frequency_mhz[indices_array]
            spacings = np.diff(guides)
            if np.all(
                (spacings >= minimum_spacing)
                & (spacings <= maximum_spacing)
            ):
                chains.append((indices_array, guides))
    return chains


def track_path_frequency_score(
    waterfall: np.ndarray,
    frequency_mhz: np.ndarray,
    guides_mhz: np.ndarray,
    result: dict[str, np.ndarray],
) -> tuple[float, float, float]:
    guide_channels = nearest_channels(frequency_mhz, guides_mhz)
    time_count = waterfall.shape[0]
    rows = np.arange(time_count)[:, None, None]
    local_offsets = np.arange(-2, 3, dtype=np.int64)[None, None, :]

    def score(common_offset: int) -> float:
        channels = (
            guide_channels[None, :, None]
            + result["common_shift_channels"][:, None, None]
            + common_offset
            + local_offsets
        )
        if np.min(channels) < 0 or np.max(channels) >= waterfall.shape[1]:
            return np.nan
        local_peak = np.nanmax(waterfall[rows, channels], axis=2)
        return float(np.nanmedian(np.nanmedian(local_peak, axis=1)))

    fixed_score = score(0)
    null_offsets = list(range(-18, -6)) + list(range(7, 19))
    null_scores = np.asarray([score(offset) for offset in null_offsets])
    null_center = float(np.nanmedian(null_scores))
    null_mad = float(np.nanmedian(np.abs(null_scores - null_center)))
    score_z = (
        (fixed_score - null_center) / (1.4826 * null_mad)
        if null_mad > 0
        else np.nan
    )
    percentile = float(100.0 * np.nanmean(null_scores < fixed_score))
    return fixed_score, score_z, percentile


def compare_with_reference(
    target_mjd: np.ndarray,
    target_shift_mhz: np.ndarray,
    reference_mjd: np.ndarray,
    reference_shift_mhz: np.ndarray,
) -> tuple[float, float, float, np.ndarray, np.ndarray]:
    overlap = (
        (target_mjd >= np.nanmin(reference_mjd))
        & (target_mjd <= np.nanmax(reference_mjd))
    )
    reference = np.interp(
        target_mjd[overlap], reference_mjd, reference_shift_mhz
    )
    target = target_shift_mhz[overlap]
    difference_khz = (target - reference) * 1000.0
    median_offset_khz = float(np.nanmedian(difference_khz))
    residual_mad_khz = float(
        np.nanmedian(np.abs(difference_khz - median_offset_khz))
    )
    return (
        correlation(target, reference),
        median_offset_khz,
        residual_mad_khz,
        overlap,
        reference,
    )


def fit_reference_assisted_track(
    waterfall: np.ndarray,
    frequency_mhz: np.ndarray,
    target_mjd: np.ndarray,
    guides_mhz: np.ndarray,
    reference_mjd: np.ndarray,
    reference_shift_mhz: np.ndarray,
    support_threshold: float,
) -> dict[str, np.ndarray | float]:
    overlap = (
        (target_mjd >= np.nanmin(reference_mjd))
        & (target_mjd <= np.nanmax(reference_mjd))
    )
    target_rows = np.flatnonzero(overlap)
    reference = np.interp(
        target_mjd[overlap], reference_mjd, reference_shift_mhz
    )
    delta_frequency_mhz = float(np.median(np.diff(frequency_mhz)))
    reference_channels = np.rint(
        reference / delta_frequency_mhz
    ).astype(np.int64)
    guide_channels = nearest_channels(frequency_mhz, guides_mhz)
    common_offsets = np.arange(-18, 19, dtype=np.int64)
    local_offsets = np.arange(-2, 3, dtype=np.int64)
    response = waterfall[overlap]
    time_indices = np.arange(len(response))[:, None, None]
    scores = np.full(
        (len(response), len(common_offsets)), np.nan, dtype=np.float64
    )
    supports = np.zeros_like(scores, dtype=np.int64)

    for state_index, common_offset in enumerate(common_offsets):
        channels = (
            guide_channels[None, :, None]
            + reference_channels[:, None, None]
            + common_offset
            + local_offsets[None, None, :]
        )
        if np.min(channels) < 0 or np.max(channels) >= response.shape[1]:
            continue
        local_peak = np.nanmax(response[time_indices, channels], axis=2)
        clipped = np.sort(np.clip(local_peak, -2.0, 8.0), axis=1)
        used_teeth = min(3, len(guides_mhz))
        scores[:, state_index] = np.nanmean(
            clipped[:, -used_teeth:], axis=1
        )
        supports[:, state_index] = np.count_nonzero(
            local_peak >= support_threshold, axis=1
        )

    aggregate = np.nanmedian(scores, axis=0)
    best_state = int(np.nanargmax(aggregate))
    selected_score = scores[:, best_state]
    row_center = np.nanmedian(scores, axis=1)
    row_scale = robust_scale(scores, axis=1)
    confidence = (selected_score - row_center) / row_scale
    minimum_support = min(3, len(guides_mhz))
    reliable = (
        (supports[:, best_state] >= minimum_support)
        & (confidence >= 1.0)
    )
    aggregate_center = float(np.nanmedian(aggregate))
    aggregate_scale = float(robust_scale(aggregate))
    aggregate_z = float(
        (aggregate[best_state] - aggregate_center) / aggregate_scale
    )
    fitted_offset_channels = int(common_offsets[best_state])
    fitted_shift_mhz = (
        reference
        + fitted_offset_channels * delta_frequency_mhz
    )
    return {
        "overlap_rows": target_rows,
        "shift_mhz": fitted_shift_mhz,
        "track_mhz": guides_mhz[None, :] + fitted_shift_mhz[:, None],
        "confidence": confidence,
        "support_count": supports[:, best_state],
        "reliable": reliable,
        "fitted_offset_channels": fitted_offset_channels,
        "fitted_offset_khz": (
            fitted_offset_channels * delta_frequency_mhz * 1000.0
        ),
        "aggregate_score_z": aggregate_z,
    }


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_best_result(
    path: Path,
    scan: int,
    beam: str,
    polar: str,
    frequency_mhz: np.ndarray,
    elapsed_minutes: np.ndarray,
    profile: np.ndarray,
    waterfall: np.ndarray,
    guides_mhz: np.ndarray,
    result: dict[str, np.ndarray],
    reference_assisted: dict[str, np.ndarray | float],
    accepted: bool,
) -> None:
    figure, axes = plt.subplots(
        3,
        1,
        figsize=(14, 10),
        gridspec_kw={"height_ratios": [2.0, 0.7, 0.8]},
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
        overlap_rows = reference_assisted["overlap_rows"]
        axes[0].plot(
            reference_assisted["track_mhz"][:, tooth_index],
            elapsed_minutes[overlap_rows],
            color="magenta",
            linewidth=0.9,
            linestyle="--",
        )
    axes[0].set(
        ylabel="time from start (min)",
        title=(
            f"scan {scan} {beam}/{polar} local guides; "
            "black=independent, magenta=reference-assisted"
        ),
    )
    figure.colorbar(
        image,
        ax=axes[0],
        pad=0.01,
        label="smoothed response / robust RMS",
    )

    axes[1].plot(frequency_mhz, profile, color="tab:blue")
    axes[1].scatter(
        guides_mhz,
        np.interp(guides_mhz, frequency_mhz, profile),
        color="black",
        zorder=3,
    )
    axes[1].set(
        ylabel="time-summary response",
        title=f"selected local guides; accepted={accepted}",
    )
    axes[1].grid(alpha=0.2)

    axes[2].plot(
        elapsed_minutes,
        result["common_shift_mhz"] * 1000.0,
        color="black",
        label="independent tracker",
    )
    overlap_rows = reference_assisted["overlap_rows"]
    axes[2].plot(
        elapsed_minutes[overlap_rows],
        reference_assisted["shift_mhz"] * 1000.0,
        color="magenta",
        linestyle="--",
        label="reference-assisted",
    )
    axes[2].set(
        xlabel="time from start (min)",
        ylabel="common shift (kHz)",
    )
    axes[2].grid(alpha=0.2)
    axes[2].legend()

    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> int:
    arguments = parse_arguments()
    reference_polar = arguments.reference_polar or arguments.polar
    with arguments.config.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    discovery_config = config["local_guide_discovery"]
    ranking_config = config["candidate_ranking"]
    acceptance_config = config["acceptance"]
    base_parameters = load_tracker_parameters(config["tracker"])

    if arguments.independent_only:
        reference_mjd = None
        reference_shift_mhz = None
    else:
        reference_track_path = arguments.tables_dir / (
            f"dynamic_mask_scan{arguments.scan}_"
            f"{arguments.reference_beam}_{reference_polar}_"
            f"{arguments.reference_variant}_tracks.npz"
        )
        with np.load(reference_track_path) as product:
            if "binned_mjd" not in product:
                reference_waterfall_path = arguments.tables_dir / (
                    f"dynamic_mask_scan{arguments.scan}_"
                    f"{arguments.reference_beam}_{reference_polar}_"
                    f"{arguments.reference_variant}_waterfall.npz"
                )
                with np.load(
                    reference_waterfall_path
                ) as waterfall_product:
                    reference_mjd = waterfall_product[
                        "binned_mjd"
                    ].astype(np.float64)
            else:
                reference_mjd = product["binned_mjd"].astype(np.float64)
            shift_key = (
                "common_shift_mhz"
                if "common_shift_mhz" in product
                else "independent_common_shift_mhz"
            )
            reference_shift_mhz = product[shift_key].astype(np.float64)

    all_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for beam in arguments.beams:
        waterfall_path = arguments.tables_dir / (
            f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
            f"{arguments.waterfall_variant}_waterfall.npz"
        )
        with np.load(waterfall_path) as product:
            frequency_mhz = product["frequency_mhz"].astype(np.float64)
            waterfall = product["smoothed_standardized"].astype(np.float64)
            target_mjd = product["binned_mjd"].astype(np.float64)
            elapsed_minutes = product[
                "binned_elapsed_minutes"
            ].astype(np.float64)

        profile = build_profile(waterfall)
        peaks = profile_peaks(frequency_mhz, profile, discovery_config)
        if arguments.exclude_mask_product is not None:
            with np.load(
                arguments.exclude_mask_product, allow_pickle=False
            ) as exclusion_product:
                excluded_centers = exclusion_product[
                    "fitted_centers_mhz"
                ].astype(np.float64)
            keep = np.ones(len(peaks), dtype=bool)
            for center_mhz in excluded_centers:
                keep &= (
                    np.abs(frequency_mhz[peaks] - center_mhz)
                    > arguments.exclude_half_width_mhz
                )
            peaks = peaks[keep]
        chains = candidate_chains(
            frequency_mhz, peaks, discovery_config
        )
        evaluated = []
        for profile_indices, guides_mhz in chains:
            tooth_count = len(guides_mhz)
            parameters = replace(
                base_parameters,
                guide_centers_mhz=guides_mhz,
                teeth_used_for_score=min(
                    base_parameters.teeth_used_for_score, tooth_count
                ),
                minimum_supported_teeth=min(
                    base_parameters.minimum_supported_teeth, tooth_count
                ),
            )
            result = track_comb(
                frequency_mhz, waterfall, parameters
            )
            (
                fixed_path_score,
                dynamic_path_z,
                dynamic_path_percentile,
            ) = track_path_frequency_score(
                waterfall, frequency_mhz, guides_mhz, result
            )
            if arguments.independent_only:
                reference_correlation = 0.0
                reference_offset_khz = np.nan
                reference_residual_mad_khz = 0.0
            else:
                (
                    reference_correlation,
                    reference_offset_khz,
                    reference_residual_mad_khz,
                    _,
                    _,
                ) = compare_with_reference(
                    target_mjd,
                    result["common_shift_mhz"],
                    reference_mjd,
                    reference_shift_mhz,
                )
            maximum_state = int(
                np.max(np.abs(result["state_offsets"]))
            )
            search_edge_fraction = float(
                np.mean(
                    np.abs(result["common_shift_channels"])
                    >= maximum_state
                )
            )
            profile_strength = float(
                np.nanmedian(profile[profile_indices])
            )
            length_bonus = 1 if tooth_count == 4 else 0
            clipped_dynamic_z = np.clip(
                dynamic_path_z,
                -float(ranking_config["dynamic_path_z_clip"]),
                float(ranking_config["dynamic_path_z_clip"]),
            )
            selection_score = float(
                clipped_dynamic_z
                + float(
                    ranking_config["reference_correlation_weight"]
                )
                * reference_correlation
                - float(
                    ranking_config[
                        "reference_residual_mad_khz_penalty"
                    ]
                )
                * reference_residual_mad_khz
                + float(ranking_config["four_tooth_bonus"])
                * length_bonus
                + float(ranking_config["profile_strength_weight"])
                * profile_strength
                - float(
                    ranking_config["search_edge_fraction_penalty"]
                )
                * search_edge_fraction
            )
            accepted = bool(
                dynamic_path_z
                >= float(
                    acceptance_config["minimum_dynamic_path_z"]
                )
                and (
                    arguments.independent_only
                    or (
                        reference_correlation
                        >= float(
                            acceptance_config[
                                "minimum_reference_track_correlation"
                            ]
                        )
                        and reference_residual_mad_khz
                        <= float(
                            acceptance_config[
                                "maximum_reference_track_residual_mad_khz"
                            ]
                        )
                    )
                )
                and search_edge_fraction
                <= float(
                    acceptance_config[
                        "maximum_search_edge_fraction"
                    ]
                )
            )
            evaluated.append(
                {
                    "beam": beam,
                    "tooth_count": tooth_count,
                    "guide_centers_mhz": ";".join(
                        f"{value:.9f}" for value in guides_mhz
                    ),
                    "spacings_mhz": ";".join(
                        f"{value:.9f}" for value in np.diff(guides_mhz)
                    ),
                    "profile_strength": profile_strength,
                    "fixed_path_score": fixed_path_score,
                    "dynamic_path_z": dynamic_path_z,
                    "dynamic_path_percentile": dynamic_path_percentile,
                    "reliable_time_fraction": float(
                        np.mean(result["reliable"])
                    ),
                    "search_edge_fraction": search_edge_fraction,
                    "reference_track_correlation": (
                        reference_correlation
                    ),
                    "reference_track_median_offset_khz": (
                        reference_offset_khz
                    ),
                    "reference_track_residual_mad_khz": (
                        reference_residual_mad_khz
                    ),
                    "selection_score": selection_score,
                    "accepted": accepted,
                    "_guides": guides_mhz,
                    "_result": result,
                }
            )

        if not evaluated:
            raise RuntimeError(f"No local guide chains found for {beam}")
        evaluated.sort(
            key=lambda row: float(row["selection_score"]), reverse=True
        )
        for rank, row in enumerate(evaluated, start=1):
            output_row = {
                "beam": row["beam"],
                "rank": rank,
                **{
                    key: value
                    for key, value in row.items()
                    if not key.startswith("_") and key != "beam"
                },
            }
            all_rows.append(output_row)

        if not 1 <= arguments.output_candidate_rank <= len(evaluated):
            raise ValueError(
                f"Requested output candidate rank "
                f"{arguments.output_candidate_rank} is unavailable for {beam}"
            )
        best = evaluated[arguments.output_candidate_rank - 1]
        best_guides = best["_guides"]
        best_result = best["_result"]
        if arguments.independent_only:
            reference_assisted = {
                "overlap_rows": np.asarray([], dtype=np.int64),
                "shift_mhz": np.asarray([], dtype=np.float64),
                "track_mhz": np.empty(
                    (0, len(best_guides)), dtype=np.float64
                ),
                "confidence": np.asarray([], dtype=np.float64),
                "support_count": np.asarray([], dtype=np.int64),
                "reliable": np.asarray([], dtype=bool),
                "fitted_offset_channels": 0,
                "fitted_offset_khz": np.nan,
                "aggregate_score_z": np.nan,
            }
        else:
            reference_assisted = fit_reference_assisted_track(
                waterfall,
                frequency_mhz,
                target_mjd,
                best_guides,
                reference_mjd,
                reference_shift_mhz,
                base_parameters.support_response_threshold,
            )
        summary_rows.append(
            {
                "beam": beam,
                "rank": arguments.output_candidate_rank,
                **{
                    key: value
                    for key, value in best.items()
                    if not key.startswith("_") and key != "beam"
                },
            }
            | {
                "reference_assisted_offset_khz": (
                    reference_assisted["fitted_offset_khz"]
                ),
                "reference_assisted_score_z": (
                    reference_assisted["aggregate_score_z"]
                ),
                "reference_assisted_reliable_fraction": float(
                    np.mean(reference_assisted["reliable"])
                    if len(reference_assisted["reliable"])
                    else np.nan
                ),
            }
        )

        output_stem = (
            f"dynamic_mask_scan{arguments.scan}_{beam}_{arguments.polar}_"
            f"{arguments.output_label}"
        )
        output_npz = arguments.tables_dir / f"{output_stem}_tracks.npz"
        np.savez_compressed(
            output_npz,
            beam=beam,
            accepted=best["accepted"],
            frequency_mhz=frequency_mhz,
            binned_mjd=target_mjd,
            binned_elapsed_minutes=elapsed_minutes,
            guide_centers_mhz=best_guides,
            profile=profile,
            independent_common_shift_mhz=best_result[
                "common_shift_mhz"
            ],
            independent_common_track_mhz=best_result[
                "common_track_mhz"
            ],
            independent_reliable=best_result["reliable"],
            reference_assisted_overlap_rows=reference_assisted[
                "overlap_rows"
            ],
            reference_assisted_shift_mhz=reference_assisted[
                "shift_mhz"
            ],
            reference_assisted_track_mhz=reference_assisted[
                "track_mhz"
            ],
            reference_assisted_reliable=reference_assisted[
                "reliable"
            ],
            reference_assisted_confidence=reference_assisted[
                "confidence"
            ],
            reference_assisted_support_count=reference_assisted[
                "support_count"
            ],
        )
        plot_best_result(
            arguments.figures_dir / f"{output_stem}.png",
            arguments.scan,
            beam,
            arguments.polar,
            frequency_mhz,
            elapsed_minutes,
            profile,
            waterfall,
            best_guides,
            best_result,
            reference_assisted,
            bool(best["accepted"]),
        )

    write_rows(arguments.candidate_csv, all_rows)
    write_rows(arguments.summary_csv, summary_rows)
    print(f"saved_candidate_csv={arguments.candidate_csv}")
    print(f"saved_summary_csv={arguments.summary_csv}")
    for row in summary_rows:
        print(
            f"{row['beam']}: accepted={row['accepted']}, "
            f"guides={row['guide_centers_mhz']}, "
            f"reliable={float(row['reliable_time_fraction']):.6f}, "
            f"reference_r={float(row['reference_track_correlation']):.6f}, "
            "reference_residual_mad_khz="
            f"{float(row['reference_track_residual_mad_khz']):.6f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
