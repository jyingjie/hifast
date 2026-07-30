#!/usr/bin/env python3
"""Blindly detect local positive-bump sequences in group-x spectra."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import find_peaks, peak_widths, savgol_filter


WIDTH_BANK_MHZ = [0.06, 0.08, 0.10, 0.14, 0.18, 0.22, 0.26]


@dataclass(frozen=True)
class Bump:
    channel: int
    frequency_mhz: float
    fwhm_mhz: float
    matched_fwhm_mhz: float
    standardized_response: float
    prominence_jy_beam: float


def robust_scale(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if not len(finite):
        return float("nan")
    center = np.median(finite)
    return float(1.4826 * np.median(np.abs(finite - center)))


def odd_window(width_mhz: float, delta_freq_mhz: float) -> int:
    value = max(3, int(round(width_mhz / abs(delta_freq_mhz))))
    return value if value % 2 else value + 1


def fill_spectrum(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    if np.count_nonzero(finite) < 2:
        return np.full_like(values, np.nan)
    result = values.copy()
    indices = np.arange(len(values))
    result[~finite] = np.interp(
        indices[~finite], indices[finite], values[finite]
    )
    return result


def baseline_residual(
    values: np.ndarray,
    delta_freq_mhz: float,
    method: str,
) -> np.ndarray:
    filled = fill_spectrum(values)
    if not np.all(np.isfinite(filled)):
        return np.full_like(values, np.nan)
    if method == "local_median":
        baseline = median_filter(
            filled,
            size=odd_window(0.55, delta_freq_mhz),
            mode="nearest",
        )
    elif method == "savgol":
        baseline = savgol_filter(
            filled,
            window_length=odd_window(0.65, delta_freq_mhz),
            polyorder=2,
            mode="interp",
        )
    else:
        raise ValueError(method)
    return values - baseline


def detect_bumps(
    freq: np.ndarray,
    residual: np.ndarray,
    *,
    sign: int = 1,
    min_response: float = 1.5,
    min_fwhm_mhz: float = 0.06,
    max_fwhm_mhz: float = 0.26,
) -> list[Bump]:
    delta_freq = float(np.median(np.diff(freq)))
    values = sign * fill_spectrum(residual)
    if not np.all(np.isfinite(values)):
        return []

    standardized = []
    for width in WIDTH_BANK_MHZ:
        sigma_chan = width / 2.354820045 / abs(delta_freq)
        response = gaussian_filter1d(values, sigma=sigma_chan, mode="nearest")
        scale = robust_scale(response)
        if not np.isfinite(scale) or scale <= 0:
            standardized.append(np.full_like(response, np.nan))
        else:
            standardized.append(response / scale)
    response_stack = np.stack(standardized)
    best_response = np.nanmax(response_stack, axis=0)
    best_width_index = np.nanargmax(response_stack, axis=0)

    matched_peaks, _ = find_peaks(
        best_response,
        distance=max(1, int(round(0.11 / abs(delta_freq)))),
        prominence=0.4,
    )

    lightly_smoothed = gaussian_filter1d(values, sigma=1.5, mode="nearest")
    shape_peaks, shape_properties = find_peaks(
        lightly_smoothed, prominence=0
    )
    if not len(shape_peaks):
        return []
    shape_widths = (
        peak_widths(
            lightly_smoothed, shape_peaks, rel_height=0.5
        )[0]
        * abs(delta_freq)
    )

    edge_mhz = 0.30
    by_shape_channel: dict[int, Bump] = {}
    for matched_channel in matched_peaks:
        nearest_index = int(
            np.argmin(np.abs(shape_peaks - matched_channel))
        )
        shape_channel = int(shape_peaks[nearest_index])
        if abs(shape_channel - matched_channel) * abs(delta_freq) > 0.08:
            continue
        measured_width = float(shape_widths[nearest_index])
        score = float(best_response[matched_channel])
        frequency = float(freq[shape_channel])
        if not (
            min_response <= score
            and min_fwhm_mhz <= measured_width <= max_fwhm_mhz
            and freq[0] + edge_mhz <= frequency <= freq[-1] - edge_mhz
        ):
            continue
        bump = Bump(
            channel=shape_channel,
            frequency_mhz=frequency,
            fwhm_mhz=measured_width,
            matched_fwhm_mhz=WIDTH_BANK_MHZ[
                int(best_width_index[matched_channel])
            ],
            standardized_response=score,
            prominence_jy_beam=float(
                shape_properties["prominences"][nearest_index]
            ),
        )
        previous = by_shape_channel.get(shape_channel)
        if (
            previous is None
            or bump.standardized_response > previous.standardized_response
        ):
            by_shape_channel[shape_channel] = bump
    return sorted(by_shape_channel.values(), key=lambda bump: bump.frequency_mhz)


def best_sequence(
    bumps: list[Bump],
    minimum_spacing_mhz: float = 0.80,
    maximum_spacing_mhz: float = 1.10,
) -> list[Bump]:
    if not bumps:
        return []
    paths: list[list[int]] = [[index] for index in range(len(bumps))]
    for right in range(len(bumps)):
        for left in range(right):
            spacing = (
                bumps[right].frequency_mhz - bumps[left].frequency_mhz
            )
            if minimum_spacing_mhz <= spacing <= maximum_spacing_mhz:
                proposal = paths[left] + [right]
                current = paths[right]
                proposal_score = sum(
                    bumps[index].standardized_response
                    for index in proposal
                )
                current_score = sum(
                    bumps[index].standardized_response
                    for index in current
                )
                if (
                    len(proposal) > len(current)
                    or (
                        len(proposal) == len(current)
                        and proposal_score > current_score
                    )
                ):
                    paths[right] = proposal
    best_path = max(
        paths,
        key=lambda path: (
            len(path),
            sum(bumps[index].standardized_response for index in path),
        ),
    )
    return [bumps[index] for index in best_path]


def confirm_chain_with_alternate_baseline(
    freq: np.ndarray,
    residual: np.ndarray,
    reference_chain: list[Bump],
    minimum_response: float = 1.5,
) -> list[Bump]:
    """Measure a previously blind-detected chain after changing the baseline."""
    delta_freq = float(np.median(np.diff(freq)))
    confirmations = []
    for reference in reference_chain:
        search = np.flatnonzero(
            np.abs(freq - reference.frequency_mhz) <= 0.08
        )
        best = None
        for width in WIDTH_BANK_MHZ:
            smoothed = gaussian_filter1d(
                residual,
                sigma=width / 2.354820045 / abs(delta_freq),
                mode="nearest",
            )
            noise = robust_scale(smoothed)
            channel = int(search[np.argmax(smoothed[search] / noise)])
            score = float(smoothed[channel] / noise)
            if best is None or score > best.standardized_response:
                best = Bump(
                    channel=channel,
                    frequency_mhz=float(freq[channel]),
                    fwhm_mhz=float(width),
                    matched_fwhm_mhz=float(width),
                    standardized_response=score,
                    prominence_jy_beam=float("nan"),
                )
        if best is not None and best.standardized_response >= minimum_response:
            confirmations.append(best)
    return confirmations


def aggregate_blocks(
    spectra: np.ndarray,
    weights: np.ndarray,
    factor: int,
) -> tuple[np.ndarray, np.ndarray]:
    group_count, _, x_count, channel_count = spectra.shape
    block_count = (x_count + factor - 1) // factor
    result = np.full(
        (group_count, block_count, channel_count),
        np.nan,
        dtype=np.float64,
    )
    block_bounds = np.zeros((block_count, 2), dtype=int)
    merged = 0.5 * (spectra[:, 0] + spectra[:, 1])
    merged_weight = np.minimum(weights[:, 0], weights[:, 1])
    for block in range(block_count):
        start = block * factor
        stop = min((block + 1) * factor, x_count)
        block_bounds[block] = start, stop
        values = merged[:, start:stop].astype(np.float64)
        use_weight = merged_weight[:, start:stop].astype(np.float64)
        numerator = np.nansum(values * use_weight, axis=1)
        denominator = np.nansum(
            np.where(np.isfinite(values), use_weight, 0),
            axis=1,
        )
        result[:, block] = np.divide(
            numerator,
            denominator,
            out=np.full_like(numerator, np.nan),
            where=denominator > 0,
        )
    return result, block_bounds


def chain_summary(chain: list[Bump]) -> tuple[str, str, str, str, float]:
    centers = ";".join(f"{bump.frequency_mhz:.9f}" for bump in chain)
    widths = ";".join(f"{bump.fwhm_mhz:.9f}" for bump in chain)
    responses = ";".join(
        f"{bump.standardized_response:.9g}" for bump in chain
    )
    spacings_values = [
        right.frequency_mhz - left.frequency_mhz
        for left, right in zip(chain[:-1], chain[1:])
    ]
    spacings = ";".join(f"{value:.9f}" for value in spacings_values)
    score = float(sum(bump.standardized_response for bump in chain))
    return centers, widths, responses, spacings, score


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def permute_chunks(
    values: np.ndarray,
    rng: np.random.Generator,
    chunk_channels: int,
) -> np.ndarray:
    chunks = [
        values[start : min(start + chunk_channels, len(values))]
        for start in range(0, len(values), chunk_channels)
    ]
    order = rng.permutation(len(chunks))
    return np.concatenate([chunks[index] for index in order])


def time_shift_frequency_chunks(
    values: np.ndarray,
    weights: np.ndarray,
    rng: np.random.Generator,
    chunk_channels: int,
    block_start: int,
    block_stop: int,
) -> np.ndarray:
    """Break cross-frequency time coincidence using independent x shifts."""
    numerator_parts = []
    denominator_parts = []
    for channel_start in range(0, values.shape[1], chunk_channels):
        channel_stop = min(channel_start + chunk_channels, values.shape[1])
        shift = int(rng.integers(values.shape[0]))
        shifted_values = np.roll(
            values[:, channel_start:channel_stop], shift, axis=0
        )[block_start:block_stop]
        shifted_weights = np.roll(
            weights[:, channel_start:channel_stop], shift, axis=0
        )[block_start:block_stop]
        numerator_parts.append(
            np.nansum(shifted_values * shifted_weights, axis=0)
        )
        denominator_parts.append(
            np.nansum(
                np.where(np.isfinite(shifted_values), shifted_weights, 0.0),
                axis=0,
            )
        )
    numerator = np.concatenate(numerator_parts)
    denominator = np.concatenate(denominator_parts)
    return np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=denominator > 0,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--block-factor", type=int, default=18)
    parser.add_argument("--affected-scan", type=int, required=True)
    parser.add_argument("--affected-date", required=True)
    parser.add_argument("--affected-source", default="main")
    parser.add_argument("--affected-x-range", type=float, nargs=2, required=True)
    parser.add_argument("--control-scan", type=int, required=True)
    parser.add_argument("--control-date", required=True)
    parser.add_argument("--control-source", default="main")
    parser.add_argument("--shuffle-count", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260726)
    parser.add_argument("--bumps-csv", type=Path, required=True)
    parser.add_argument("--sequences-csv", type=Path, required=True)
    parser.add_argument("--null-tests-csv", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()

    product = np.load(args.input)
    freq = product["freq_lsrk_mhz"]
    x_centers = product["x_centers"]
    spectra = product["group_mean_spectra_jy_beam"]
    weights = product["group_finite_weight"]
    scans = product["group_scan"]
    dates = product["group_date"]
    sources = product["group_source"]
    delta_freq = float(np.median(np.diff(freq)))

    block_spectra, block_bounds = aggregate_blocks(
        spectra, weights, args.block_factor
    )
    sequence_rows = []
    bump_rows = []
    block_results = {}
    for group in range(len(scans)):
        for block in range(block_spectra.shape[1]):
            values = block_spectra[group, block]
            residual = baseline_residual(
                values, delta_freq, "local_median"
            )
            bumps = detect_bumps(freq, residual)
            chain = best_sequence(bumps)
            start, stop = block_bounds[block]
            x_min = float(x_centers[start])
            x_max = float(x_centers[stop - 1])
            key = (group, block)
            block_results[key] = (values, residual, bumps, chain)
            for candidate_number, bump in enumerate(bumps, start=1):
                bump_rows.append(
                    {
                        "scan": int(scans[group]),
                        "date": str(dates[group]),
                        "source": str(sources[group]),
                        "block": block,
                        "x_min": f"{x_min:.6f}",
                        "x_max": f"{x_max:.6f}",
                        "candidate_number": candidate_number,
                        "frequency_mhz": f"{bump.frequency_mhz:.9f}",
                        "measured_fwhm_mhz": f"{bump.fwhm_mhz:.9f}",
                        "matched_fwhm_mhz": f"{bump.matched_fwhm_mhz:.6f}",
                        "standardized_response": (
                            f"{bump.standardized_response:.9g}"
                        ),
                        "prominence_jy_beam": (
                            f"{bump.prominence_jy_beam:.12g}"
                        ),
                    }
                )
            if len(chain) >= 2:
                centers, widths, responses, spacings, score = (
                    chain_summary(chain)
                )
                sequence_rows.append(
                    {
                        "scan": int(scans[group]),
                        "date": str(dates[group]),
                        "source": str(sources[group]),
                        "block": block,
                        "x_min": f"{x_min:.6f}",
                        "x_max": f"{x_max:.6f}",
                        "chain_length": len(chain),
                        "chain_score": f"{score:.9g}",
                        "centers_mhz": centers,
                        "spacings_mhz": spacings,
                        "measured_fwhm_mhz": widths,
                        "standardized_responses": responses,
                    }
                )

    def locate_group(scan: int, date: str, source: str) -> int:
        use = (
            (scans == scan)
            & (dates == date)
            & (sources == source)
        )
        indices = np.flatnonzero(use)
        if len(indices) != 1:
            raise RuntimeError(
                f"Could not uniquely locate {scan}/{date}/{source}"
            )
        return int(indices[0])

    affected_group = locate_group(
        args.affected_scan, args.affected_date, args.affected_source
    )
    control_group = locate_group(
        args.control_scan, args.control_date, args.control_source
    )
    affected_block_candidates = []
    control_block_candidates = []
    for block, (start, stop) in enumerate(block_bounds):
        block_x_min = x_centers[start]
        block_x_max = x_centers[stop - 1]
        overlap = (
            block_x_max >= args.affected_x_range[0]
            and block_x_min <= args.affected_x_range[1]
        )
        if overlap:
            affected_block_candidates.append(block)
            control_block_candidates.append(block)
    if len(affected_block_candidates) != 1:
        raise RuntimeError(
            "Affected x range must select exactly one coarse block; "
            f"selected {affected_block_candidates}"
        )
    affected_block = affected_block_candidates[0]
    control_block = control_block_candidates[0]

    affected_values, affected_residual, affected_bumps, affected_chain = (
        block_results[(affected_group, affected_block)]
    )
    control_values, control_residual, control_bumps, control_chain = (
        block_results[(control_group, control_block)]
    )
    negative_bumps = detect_bumps(freq, affected_residual, sign=-1)
    negative_chain = best_sequence(negative_bumps)
    affected_savgol = baseline_residual(
        affected_values, delta_freq, "savgol"
    )
    savgol_bumps = detect_bumps(freq, affected_savgol)
    savgol_chain = best_sequence(savgol_bumps)
    savgol_confirmed_chain = confirm_chain_with_alternate_baseline(
        freq, affected_savgol, affected_chain
    )

    rng = np.random.default_rng(args.seed)
    chunk_channels = odd_window(0.30, delta_freq)
    shuffle_lengths = np.zeros(args.shuffle_count, dtype=int)
    shuffle_scores = np.zeros(args.shuffle_count)
    for iteration in range(args.shuffle_count):
        shuffled = permute_chunks(
            affected_residual, rng, chunk_channels
        )
        shuffled_bumps = detect_bumps(freq, shuffled)
        shuffled_chain = best_sequence(shuffled_bumps)
        shuffle_lengths[iteration] = len(shuffled_chain)
        shuffle_scores[iteration] = sum(
            bump.standardized_response for bump in shuffled_chain
        )

    affected_start, affected_stop = block_bounds[affected_block]
    affected_x_values = 0.5 * (
        spectra[affected_group, 0] + spectra[affected_group, 1]
    )
    affected_x_weights = np.minimum(
        weights[affected_group, 0], weights[affected_group, 1]
    )
    time_rng = np.random.default_rng(args.seed + 1)
    time_shuffle_lengths = np.zeros(args.shuffle_count, dtype=int)
    time_shuffle_scores = np.zeros(args.shuffle_count)
    for iteration in range(args.shuffle_count):
        time_shuffled_values = time_shift_frequency_chunks(
            affected_x_values,
            affected_x_weights,
            time_rng,
            chunk_channels,
            int(affected_start),
            int(affected_stop),
        )
        time_shuffled_residual = baseline_residual(
            time_shuffled_values, delta_freq, "local_median"
        )
        time_shuffled_bumps = detect_bumps(freq, time_shuffled_residual)
        time_shuffled_chain = best_sequence(time_shuffled_bumps)
        time_shuffle_lengths[iteration] = len(time_shuffled_chain)
        time_shuffle_scores[iteration] = sum(
            bump.standardized_response for bump in time_shuffled_chain
        )

    affected_centers, _, _, affected_spacings, affected_score = (
        chain_summary(affected_chain)
    )
    control_score = sum(
        bump.standardized_response for bump in control_chain
    )
    negative_score = sum(
        bump.standardized_response for bump in negative_chain
    )
    savgol_score = sum(
        bump.standardized_response for bump in savgol_chain
    )
    savgol_confirmed_score = sum(
        bump.standardized_response for bump in savgol_confirmed_chain
    )
    shuffle_extreme = (
        (shuffle_lengths > len(affected_chain))
        | (
            (shuffle_lengths == len(affected_chain))
            & (shuffle_scores >= affected_score)
        )
    )
    empirical_p = (1 + np.count_nonzero(shuffle_extreme)) / (
        args.shuffle_count + 1
    )
    time_shuffle_extreme = (
        (time_shuffle_lengths > len(affected_chain))
        | (
            (time_shuffle_lengths == len(affected_chain))
            & (time_shuffle_scores >= affected_score)
        )
    )
    time_shuffle_empirical_p = (
        1 + np.count_nonzero(time_shuffle_extreme)
    ) / (args.shuffle_count + 1)

    all_lengths = np.asarray(
        [
            len(result[3])
            for result in block_results.values()
            if np.all(np.isfinite(result[1]))
        ]
    )
    all_scores = np.asarray(
        [
            sum(bump.standardized_response for bump in result[3])
            for result in block_results.values()
            if np.all(np.isfinite(result[1]))
        ]
    )
    empirical_more_extreme = (
        (all_lengths > len(affected_chain))
        | (
            (all_lengths == len(affected_chain))
            & (all_scores >= affected_score)
        )
    )

    null_rows = [
        {
            "test": "affected_local_median",
            "chain_length": len(affected_chain),
            "chain_score": f"{affected_score:.9g}",
            "centers_mhz": affected_centers,
            "spacings_mhz": affected_spacings,
            "empirical_p": "",
            "notes": "Blind positive-bump detector",
        },
        {
            "test": "same_sample_savgol_baseline",
            "chain_length": len(savgol_chain),
            "chain_score": f"{savgol_score:.9g}",
            "centers_mhz": chain_summary(savgol_chain)[0],
            "spacings_mhz": chain_summary(savgol_chain)[3],
            "empirical_p": "",
            "notes": "Baseline-method stability check",
        },
        {
            "test": "savgol_confirmation_at_blind_chain",
            "chain_length": len(savgol_confirmed_chain),
            "chain_score": f"{savgol_confirmed_score:.9g}",
            "centers_mhz": chain_summary(savgol_confirmed_chain)[0],
            "spacings_mhz": chain_summary(savgol_confirmed_chain)[3],
            "empirical_p": "",
            "notes": (
                "Local-median chain positions tested after changing "
                "the baseline; not an independent blind search"
            ),
        },
        {
            "test": "matched_position_control",
            "chain_length": len(control_chain),
            "chain_score": f"{control_score:.9g}",
            "centers_mhz": chain_summary(control_chain)[0],
            "spacings_mhz": chain_summary(control_chain)[3],
            "empirical_p": "",
            "notes": "Different scan, same x block",
        },
        {
            "test": "negative_bumps_in_affected",
            "chain_length": len(negative_chain),
            "chain_score": f"{negative_score:.9g}",
            "centers_mhz": chain_summary(negative_chain)[0],
            "spacings_mhz": chain_summary(negative_chain)[3],
            "empirical_p": "",
            "notes": "Same detector applied to negative residual",
        },
        {
            "test": "frequency_chunk_shuffle",
            "chain_length": len(affected_chain),
            "chain_score": f"{affected_score:.9g}",
            "centers_mhz": "",
            "spacings_mhz": "",
            "empirical_p": f"{empirical_p:.9g}",
            "notes": (
                f"{args.shuffle_count} permutations; "
                f"chunk≈{chunk_channels * abs(delta_freq):.3f} MHz"
            ),
        },
        {
            "test": "time_shifted_frequency_chunks",
            "chain_length": len(affected_chain),
            "chain_score": f"{affected_score:.9g}",
            "centers_mhz": "",
            "spacings_mhz": "",
            "empirical_p": f"{time_shuffle_empirical_p:.9g}",
            "notes": (
                f"{args.shuffle_count} permutations; each "
                f"≈{chunk_channels * abs(delta_freq):.3f} MHz frequency "
                "chunk independently circular-shifted along fine x bins"
            ),
        },
        {
            "test": "all_group_x_blocks",
            "chain_length": len(affected_chain),
            "chain_score": f"{affected_score:.9g}",
            "centers_mhz": "",
            "spacings_mhz": "",
            "empirical_p": (
                f"{np.mean(empirical_more_extreme):.9g}"
            ),
            "notes": (
                f"Fraction of {len(all_lengths)} valid coarse x blocks "
                "with a more extreme chain; not a false-positive p-value"
            ),
        },
    ]

    write_csv(args.bumps_csv, bump_rows)
    write_csv(args.sequences_csv, sequence_rows)
    write_csv(args.null_tests_csv, null_rows)

    figure, axes = plt.subplots(
        3, 1, figsize=(15, 11), gridspec_kw={"height_ratios": [1, 1, 0.8]}
    )
    for axis, residual, bumps, chain, title in (
        (
            axes[0],
            affected_residual,
            affected_bumps,
            affected_chain,
            "Affected sample: blind positive-bump sequence",
        ),
        (
            axes[1],
            control_residual,
            control_bumps,
            control_chain,
            "Matched-position control: same detector",
        ),
    ):
        axis.plot(freq, residual * 1000, color="0.25", lw=0.8)
        chain_channels = {bump.channel for bump in chain}
        for bump in bumps:
            color = "tab:red" if bump.channel in chain_channels else "tab:blue"
            axis.axvspan(
                bump.frequency_mhz - 0.5 * bump.fwhm_mhz,
                bump.frequency_mhz + 0.5 * bump.fwhm_mhz,
                color=color,
                alpha=0.10,
            )
            axis.plot(
                bump.frequency_mhz,
                residual[bump.channel] * 1000,
                "o",
                color=color,
                ms=4,
            )
        axis.set_ylabel("M residual (mJy/beam)")
        axis.set_title(title)
        axis.grid(alpha=0.2)
        axis.set_xlim(freq[0], freq[-1])

    axes[2].hist(
        shuffle_scores[shuffle_lengths == len(affected_chain)],
        bins=30,
        color="0.55",
        alpha=0.8,
        label=f"shuffles with chain length={len(affected_chain)}",
    )
    axes[2].axvline(
        affected_score,
        color="tab:red",
        lw=1.2,
        label=f"affected score={affected_score:.2f}",
    )
    axes[2].set_xlabel("chain score")
    axes[2].set_ylabel("shuffle count")
    axes[2].set_title(
        f"Frequency-chunk shuffle null test: empirical p={empirical_p:.4g}"
    )
    axes[2].grid(alpha=0.2)
    axes[2].legend(fontsize=8)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=170)
    plt.close(figure)

    print(f"affected_chain_length={len(affected_chain)}")
    print(f"affected_centers={affected_centers}")
    print(f"affected_spacings={affected_spacings}")
    print(f"shuffle_empirical_p={empirical_p:.9g}")
    print(f"time_shuffle_empirical_p={time_shuffle_empirical_p:.9g}")
    print(f"saved_bumps={args.bumps_csv}")
    print(f"saved_sequences={args.sequences_csv}")
    print(f"saved_null_tests={args.null_tests_csv}")
    print(f"saved_figure={args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
