#!/usr/bin/env python3
"""Run the frozen v5.2 hierarchical localization and mask evaluation."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import traceback
from collections import Counter
from pathlib import Path

import numpy as np
import yaml


def parse_sample(value: str) -> tuple[int, str, str]:
    fields = value.split("|")
    if len(fields) != 3:
        raise argparse.ArgumentTypeError(
            "--sample must be scan|beam|polar, for example 3|M04|YY"
        )
    scan, beam, polar = fields
    beam = beam.upper()
    polar = polar.upper()
    if not beam.startswith("M"):
        beam = f"M{int(beam):02d}"
    if polar not in {"XX", "YY"}:
        raise argparse.ArgumentTypeError("polar must be XX or YY")
    return int(scan), beam, polar


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--selection-csv", type=Path, required=True)
    parser.add_argument(
        "--sample", action="append", type=parse_sample, required=True
    )
    parser.add_argument("--tables-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--reports-dir", type=Path, required=True)
    parser.add_argument("--output-label", required=True)
    parser.add_argument("--output-detail-csv", type=Path, required=True)
    parser.add_argument("--output-summary-csv", type=Path, required=True)
    parser.add_argument(
        "--reuse-waterfalls",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--keep-going",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def boolean(value: object) -> bool:
    return str(value).strip().lower() == "true"


def finite_float(value: object, default: float = np.nan) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return parsed if np.isfinite(parsed) else default


def find_selection(
    rows: list[dict[str, str]], scan: int, beam: str
) -> dict[str, str]:
    beam_number = int(beam[1:])
    matches = [
        row
        for row in rows
        if int(row["scan"]) == scan and int(row["beam"]) == beam_number
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one selection row for scan {scan}/{beam}, "
            f"found {len(matches)}"
        )
    return matches[0]


def waterfall_path(
    tables_dir: Path, scan: int, beam: str, polar: str, variant: str
) -> Path:
    return tables_dir / (
        f"dynamic_mask_scan{scan}_{beam}_{polar}_{variant}_waterfall.npz"
    )


def waterfall_is_reusable(
    path: Path,
    source_hdf5: Path,
    polar: str,
    start: int,
    stop: int,
) -> bool:
    if not path.exists():
        return False
    try:
        with np.load(path, allow_pickle=False) as product:
            indices = product["selected_original_indices"]
            return bool(
                str(product["source_hdf5"]) == str(source_hdf5)
                and str(product["polarization"]) == polar
                and len(indices) == stop - start
                and int(indices[0]) == start
                and int(indices[-1]) == stop - 1
            )
    except (KeyError, OSError, ValueError):
        return False


def run_command(command: list[str], log_path: Path) -> None:
    environment = os.environ.copy()
    environment.setdefault(
        "MPLCONFIGDIR", str(Path.cwd() / ".mplconfig")
    )
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(command) + "\n")
        handle.write(completed.stdout)
        if completed.stdout and not completed.stdout.endswith("\n"):
            handle.write("\n")
    if completed.returncode:
        tail = "\n".join(completed.stdout.splitlines()[-12:])
        raise RuntimeError(
            f"Command failed with exit code {completed.returncode}:\n{tail}"
        )


def stage_failure_reason(error: Exception) -> str:
    message = str(error)
    if "No candidate pairs meet the pre-track requirements" in message:
        return "no_candidate_pairs_meet_pretrack_requirements"
    compact = " ".join(message.split())
    return f"stage_execution_failed: {compact[:240]}"


def build_waterfalls(
    *,
    source: Path,
    start: int,
    stop: int,
    scan: int,
    beam: str,
    polar: str,
    config: dict,
    tables_dir: Path,
    figures_dir: Path,
    source_dir: Path,
    reuse: bool,
    log_path: Path,
) -> None:
    waterfall_config = config["waterfall"]
    for seconds, variant in ((15.0, "pre_fc"), (30.0, "pre_fc_30s")):
        output = waterfall_path(
            tables_dir, scan, beam, polar, variant
        )
        if reuse and waterfall_is_reusable(
            output, source, polar, start, stop
        ):
            continue
        figure = figures_dir / (
            f"dynamic_mask_scan{scan}_{beam}_{polar}_{variant}_waterfall.png"
        )
        command = [
            sys.executable,
            str(source_dir / "build_dynamic_waterfall.py"),
            "--input",
            str(source),
            "--field",
            "flux",
            "--polar",
            polar,
            "--time-index-start",
            str(start),
            "--time-index-stop",
            str(stop),
            "--freq-range",
            *[
                str(value)
                for value in waterfall_config["frequency_range_mhz"]
            ],
            "--time-bin-seconds",
            str(seconds),
            "--baseline-width-mhz",
            str(waterfall_config["baseline_width_mhz"]),
            "--smooth-fwhm-mhz",
            str(waterfall_config["smooth_fwhm_mhz"]),
            "--sample-id",
            f"scan{scan}_{beam}_{polar}",
            "--output-npz",
            str(output),
            "--figure",
            str(figure),
        ]
        run_command(command, log_path)


def selected_row(path: Path) -> dict[str, str]:
    rows = [row for row in read_csv(path) if boolean(row["selected"])]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one selected row in {path}")
    return rows[0]


def stage_localization_accepted(
    row: dict[str, str], rules: dict
) -> tuple[bool, str]:
    matched = int(row["matched_count"])
    agreement = float(row["agreement_within_4_channels_fraction"])
    correlation = finite_float(row["shift_correlation"])
    if not boolean(row["both_individually_accepted"]):
        return False, "candidate_not_individually_accepted"
    if not boolean(row["both_spacing_regular"]):
        return False, "irregular_spacing"
    if matched < int(rules["minimum_matched_guides"]):
        return False, "too_few_matched_guides"
    if agreement < float(rules["minimum_track_agreement_fraction"]):
        return False, "low_15s_30s_agreement"
    minimum_correlation = rules.get("minimum_track_correlation")
    if minimum_correlation is not None and (
        not np.isfinite(correlation)
        or correlation < float(minimum_correlation)
    ):
        return False, "low_15s_30s_correlation"
    return True, "accepted"


def candidate_has_edge_evidence(
    rows: list[dict[str, str]], rules: dict
) -> bool:
    required_teeth = int(rules["required_tooth_count"])
    for row in rows:
        if int(row["tooth_count"]) != required_teeth:
            continue
        if (
            finite_float(row["dynamic_path_z"], -np.inf)
            < float(rules["minimum_dynamic_path_z"])
        ):
            continue
        if (
            finite_float(row["profile_strength"], -np.inf)
            < float(rules["minimum_profile_strength"])
        ):
            continue
        if (
            finite_float(row["search_edge_fraction"], -np.inf)
            <= float(
                rules["minimum_search_edge_fraction_exclusive"]
            )
        ):
            continue
        return True
    return False


def edge_search_triggered(
    left_candidates: Path, right_candidates: Path, rules: dict
) -> bool:
    return candidate_has_edge_evidence(
        read_csv(left_candidates), rules
    ) and candidate_has_edge_evidence(
        read_csv(right_candidates), rules
    )


def run_localization_stage(
    *,
    stage: str,
    stage_config: Path,
    rules: dict,
    scan: int,
    beam: str,
    polar: str,
    output_label: str,
    tables_dir: Path,
    figures_dir: Path,
    source_dir: Path,
    log_path: Path,
) -> dict[str, object]:
    label = f"{output_label}_{beam.lower()}_{stage}"
    candidates: dict[str, Path] = {}
    for suffix, variant in (("15s", "pre_fc"), ("30s", "pre_fc_30s")):
        candidate_csv = tables_dir / (
            f"dynamic_mask_scan{scan}_{polar}_{label}_"
            f"{suffix}_candidates.csv"
        )
        summary_csv = tables_dir / (
            f"dynamic_mask_scan{scan}_{polar}_{label}_{suffix}_summary.csv"
        )
        command = [
            sys.executable,
            str(source_dir / "discover_local_comb_guides.py"),
            "--config",
            str(stage_config),
            "--tables-dir",
            str(tables_dir),
            "--figures-dir",
            str(figures_dir),
            "--scan",
            str(scan),
            "--polar",
            polar,
            "--independent-only",
            "--waterfall-variant",
            variant,
            "--output-label",
            f"{label}_{suffix}",
            "--beams",
            beam,
            "--candidate-csv",
            str(candidate_csv),
            "--summary-csv",
            str(summary_csv),
        ]
        run_command(command, log_path)
        candidates[suffix] = candidate_csv

    pair_csv = tables_dir / (
        f"dynamic_mask_scan{scan}_{polar}_{label}_candidate_pairs.csv"
    )
    run_command(
        [
            sys.executable,
            str(source_dir / "select_stable_candidate_pair.py"),
            "--left-candidates",
            str(candidates["15s"]),
            "--right-candidates",
            str(candidates["30s"]),
            "--maximum-match-khz",
            str(rules["maximum_candidate_match_khz"]),
            "--maximum-spacing-std-mhz",
            str(rules["maximum_spacing_standard_deviation_mhz"]),
            "--output-csv",
            str(pair_csv),
        ],
        log_path,
    )
    track_csv = tables_dir / (
        f"dynamic_mask_scan{scan}_{polar}_{label}_"
        "track_selected_candidate_pairs.csv"
    )
    run_command(
        [
            sys.executable,
            str(source_dir / "select_track_consistent_candidate_pair.py"),
            "--config",
            str(stage_config),
            "--candidate-pairs",
            str(pair_csv),
            "--tables-dir",
            str(tables_dir),
            "--scan",
            str(scan),
            "--polar",
            polar,
            "--left-waterfall-variant",
            "pre_fc",
            "--right-waterfall-variant",
            "pre_fc_30s",
            "--maximum-disagreement-channels",
            "4.0",
            "--minimum-agreement-fraction",
            str(rules["minimum_track_agreement_fraction"]),
            "--minimum-matched-guides",
            str(rules["minimum_matched_guides"]),
            "--require-both-individually-accepted",
            "--require-spacing-regular",
            "--output-csv",
            str(track_csv),
        ],
        log_path,
    )
    selected = selected_row(track_csv)
    accepted, reason = stage_localization_accepted(selected, rules)
    return {
        "stage": stage,
        "label": label,
        "accepted": accepted,
        "reason": reason,
        "selected": selected,
        "track_csv": track_csv,
        "left_candidates": candidates["15s"],
        "right_candidates": candidates["30s"],
        "stage_config": stage_config,
    }


def evaluate_mask(
    *,
    stage_result: dict[str, object],
    config: dict,
    scan: int,
    beam: str,
    polar: str,
    tables_dir: Path,
    figures_dir: Path,
    reports_dir: Path,
    source_dir: Path,
    log_path: Path,
) -> dict[str, object]:
    stage = str(stage_result["stage"])
    label = str(stage_result["label"])
    stage_config = Path(stage_result["stage_config"])
    track_csv = Path(stage_result["track_csv"])
    left_label = f"{label}_selected_15s"
    right_label = f"{label}_selected_30s"
    for variant, track_label in (
        ("pre_fc", left_label),
        ("pre_fc_30s", right_label),
    ):
        run_command(
            [
                sys.executable,
                str(source_dir / "track_stable_local_guides.py"),
                "--config",
                str(stage_config),
                "--candidate-pairs",
                str(track_csv),
                "--tables-dir",
                str(tables_dir),
                "--figures-dir",
                str(figures_dir),
                "--scan",
                str(scan),
                "--polar",
                polar,
                "--waterfall-variant",
                variant,
                "--output-label",
                track_label,
                "--minimum-matched-guides",
                str(
                    config["stages"][stage]["minimum_matched_guides"]
                ),
                "--require-both-individually-accepted",
                "--require-spacing-regular",
            ],
            log_path,
        )

    consistency_csv = tables_dir / (
        f"dynamic_mask_scan{scan}_{polar}_{label}_"
        "time_bin_comparison.csv"
    )
    run_command(
        [
            sys.executable,
            str(source_dir / "compare_stable_local_tracks.py"),
            "--candidate-pairs",
            str(track_csv),
            "--tables-dir",
            str(tables_dir),
            "--scan",
            str(scan),
            "--polar",
            polar,
            "--left-label",
            left_label,
            "--right-label",
            right_label,
            "--maximum-disagreement-channels",
            "4.0",
            "--minimum-agreement-fraction",
            str(
                config["stages"][stage][
                    "minimum_track_agreement_fraction"
                ]
            ),
            "--require-spacing-regular",
            "--output-csv",
            str(consistency_csv),
        ],
        log_path,
    )

    evaluation_label = f"{label}_mask_eval"
    fixed_summary = tables_dir / (
        f"dynamic_mask_scan{scan}_{polar}_{label}_fixed_mask_summary.csv"
    )
    run_command(
        [
            sys.executable,
            str(source_dir / "evaluate_stable_group_masks.py"),
            "--track-consistency-csv",
            str(consistency_csv),
            "--tables-dir",
            str(tables_dir),
            "--figures-dir",
            str(figures_dir),
            "--reports-dir",
            str(reports_dir),
            "--scan",
            str(scan),
            "--polar",
            polar,
            "--primary-label",
            right_label,
            "--comparison-label",
            left_label,
            "--output-label",
            evaluation_label,
            "--output-summary-csv",
            str(fixed_summary),
            "--maximum-leakage-fraction",
            str(
                config["acceptance"][
                    "maximum_total_leakage_fraction_for_final_acceptance"
                ]
            ),
            "--maximum-mask-fraction",
            str(
                config["acceptance"][
                    "maximum_test_band_mask_fraction"
                ]
            ),
            "--conservative-half-width-fwhm",
            "0.80",
            "--conservative-margin-channels",
            str(config["mask"]["margin_channels"]),
        ],
        log_path,
    )
    stem = (
        f"dynamic_mask_scan{scan}_{beam}_{polar}_{evaluation_label}"
    )
    profile_csv = tables_dir / f"{stem}_profile_fits.csv"
    profile_rows = read_csv(profile_csv)
    snr = np.asarray(
        [finite_float(row["amplitude_snr"]) for row in profile_rows]
    )
    fwhm_khz = np.asarray(
        [finite_float(row["fwhm_khz"]) for row in profile_rows]
    )
    base_snr = float(config["mask"]["minimum_amplitude_snr"])
    maximum_fwhm_khz = float(
        config["mask"].get("maximum_fwhm_khz", np.inf)
    )
    selected_shape = (
        (snr >= base_snr) & (fwhm_khz <= maximum_fwhm_khz)
    )
    base_count = int(np.count_nonzero(selected_shape))
    minimum_count = int(config["mask"]["minimum_masked_tooth_count"])
    if base_count < minimum_count:
        return {
            "status": "weak_skip",
            "reason": (
                "fewer_than_three_teeth_passing_snr_and_shape"
                if np.isfinite(maximum_fwhm_khz)
                else "fewer_than_three_teeth_above_minimum_snr"
            ),
            "amplitude_snr": snr,
            "teeth_above_minimum_snr": base_count,
        }
    if stage == "edge":
        rules = config["stages"]["edge"]
        strict_count = int(
            np.count_nonzero(snr >= float(rules["strict_amplitude_snr"]))
        )
        if strict_count < int(rules["minimum_teeth_above_strict_snr"]):
            return {
                "status": "localization_failed",
                "reason": "edge_stage_strict_snr_failed",
                "amplitude_snr": snr,
                "teeth_above_minimum_snr": base_count,
                "teeth_above_strict_snr": strict_count,
            }

    grid_csv = tables_dir / (
        f"dynamic_mask_scan{scan}_{polar}_{label}_adaptive_width_grid.csv"
    )
    selection_csv = tables_dir / (
        f"dynamic_mask_scan{scan}_{polar}_{label}_"
        "adaptive_width_selection.csv"
    )
    width_values = [
        str(value) for value in config["mask"]["half_width_fwhm_grid"]
    ]
    run_command(
        [
            sys.executable,
            str(source_dir / "scan_mask_width_grid.py"),
            "--track-consistency-csv",
            str(consistency_csv),
            "--tables-dir",
            str(tables_dir),
            "--scan",
            str(scan),
            "--polar",
            polar,
            "--beams",
            beam,
            "--waterfall-variant",
            "pre_fc_30s",
            "--fit-label",
            evaluation_label,
            "--half-width-fwhm",
            *width_values,
            "--margin-channels",
            str(config["mask"]["margin_channels"]),
            "--minimum-amplitude-snr",
            str(base_snr),
            "--maximum-fwhm-khz",
            str(maximum_fwhm_khz),
            "--minimum-masked-tooth-count",
            str(minimum_count),
            "--maximum-leakage-fraction",
            str(
                config["acceptance"][
                    "maximum_total_leakage_fraction_for_width_selection"
                ]
            ),
            "--maximum-tooth-leakage-fraction",
            str(
                config["acceptance"][
                    "maximum_single_tooth_leakage_fraction"
                ]
            ),
            "--maximum-mask-fraction",
            str(
                config["acceptance"][
                    "maximum_test_band_mask_fraction"
                ]
            ),
            "--output-csv",
            str(grid_csv),
            "--output-selection-csv",
            str(selection_csv),
        ],
        log_path,
    )
    selected = read_csv(selection_csv)[0]
    accepted = boolean(selected["adaptive_mask_accepted"])
    return {
        "status": "masked" if accepted else "mask_failed",
        "reason": "accepted" if accepted else "mask_limits_not_met",
        "amplitude_snr": snr,
        "teeth_above_minimum_snr": base_count,
        "adaptive_selection": selected,
    }


def empty_result(
    scan: int, beam: str, polar: str, date: str
) -> dict[str, object]:
    return {
        "scan": scan,
        "date": date,
        "beam": beam,
        "polar": polar,
        "status": "",
        "selected_stage": "",
        "reason": "",
        "matched_guide_count": "",
        "track_agreement_fraction": "",
        "track_correlation": "",
        "left_dynamic_path_z": "",
        "right_dynamic_path_z": "",
        "left_search_edge_fraction": "",
        "right_search_edge_fraction": "",
        "teeth_above_minimum_snr": "",
        "amplitude_snr": "",
        "half_width_fwhm": "",
        "masked_tooth_count": "",
        "masked_fraction": "",
        "estimated_leakage_fraction": "",
        "maximum_tooth_leakage_fraction": "",
        "edge_triggered": False,
    }


def process_sample(
    *,
    sample: tuple[int, str, str],
    selection_rows: list[dict[str, str]],
    config: dict,
    config_path: Path,
    arguments: argparse.Namespace,
    source_dir: Path,
) -> dict[str, object]:
    scan, beam, polar = sample
    selection = find_selection(selection_rows, scan, beam)
    result = empty_result(scan, beam, polar, selection["date"])
    source = Path(selection["stage_path"])
    start = int(selection["time_index_start"])
    stop = int(selection["time_index_stop_exclusive"])
    sample_label = f"scan{scan}_{beam}_{polar}"
    log_path = arguments.reports_dir / (
        f"{arguments.output_label}_{sample_label}.log"
    )
    if log_path.exists():
        log_path.unlink()
    build_waterfalls(
        source=source,
        start=start,
        stop=stop,
        scan=scan,
        beam=beam,
        polar=polar,
        config=config,
        tables_dir=arguments.tables_dir,
        figures_dir=arguments.figures_dir,
        source_dir=source_dir,
        reuse=arguments.reuse_waterfalls,
        log_path=log_path,
    )

    last_stage: dict[str, object] | None = None
    expanded_result: dict[str, object] | None = None
    for stage in ("primary", "expanded"):
        stage_config = (
            config_path.parent / config["stage_configs"][stage]
        ).resolve()
        try:
            current = run_localization_stage(
                stage=stage,
                stage_config=stage_config,
                rules=config["stages"][stage],
                scan=scan,
                beam=beam,
                polar=polar,
                output_label=arguments.output_label,
                tables_dir=arguments.tables_dir,
                figures_dir=arguments.figures_dir,
                source_dir=source_dir,
                log_path=log_path,
            )
        except Exception as error:
            failed_label = (
                f"{arguments.output_label}_{beam.lower()}_{stage}"
            )
            current = {
                "stage": stage,
                "accepted": False,
                "reason": stage_failure_reason(error),
                "left_candidates": arguments.tables_dir
                / (
                    f"dynamic_mask_scan{scan}_{polar}_{failed_label}_"
                    "15s_candidates.csv"
                ),
                "right_candidates": arguments.tables_dir
                / (
                    f"dynamic_mask_scan{scan}_{polar}_{failed_label}_"
                    "30s_candidates.csv"
                ),
            }
        last_stage = current
        if stage == "expanded":
            expanded_result = current
        if current["accepted"]:
            break

    if last_stage is not None and not last_stage["accepted"]:
        if (
            expanded_result is not None
            and "left_candidates" in expanded_result
            and Path(expanded_result["left_candidates"]).exists()
            and Path(expanded_result["right_candidates"]).exists()
            and edge_search_triggered(
                Path(expanded_result["left_candidates"]),
                Path(expanded_result["right_candidates"]),
                config["edge_trigger"],
            )
        ):
            result["edge_triggered"] = True
            stage = "edge"
            stage_config = (
                config_path.parent / config["stage_configs"][stage]
            ).resolve()
            try:
                last_stage = run_localization_stage(
                    stage=stage,
                    stage_config=stage_config,
                    rules=config["stages"][stage],
                    scan=scan,
                    beam=beam,
                    polar=polar,
                    output_label=arguments.output_label,
                    tables_dir=arguments.tables_dir,
                    figures_dir=arguments.figures_dir,
                    source_dir=source_dir,
                    log_path=log_path,
                )
            except Exception as error:
                last_stage = {
                    "stage": stage,
                    "accepted": False,
                    "reason": stage_failure_reason(error),
                }

    if last_stage is None or not last_stage["accepted"]:
        result["status"] = "localization_failed"
        result["selected_stage"] = (
            str(last_stage["stage"]) if last_stage else ""
        )
        result["reason"] = (
            str(last_stage["reason"]) if last_stage else "no_stage_ran"
        )
        if last_stage is not None and "selected" in last_stage:
            selected = last_stage["selected"]
            result.update(
                {
                    "matched_guide_count": selected["matched_count"],
                    "track_agreement_fraction": selected[
                        "agreement_within_4_channels_fraction"
                    ],
                    "track_correlation": selected["shift_correlation"],
                    "left_dynamic_path_z": selected[
                        "left_dynamic_path_z"
                    ],
                    "right_dynamic_path_z": selected[
                        "right_dynamic_path_z"
                    ],
                    "left_search_edge_fraction": selected[
                        "left_search_edge_fraction"
                    ],
                    "right_search_edge_fraction": selected[
                        "right_search_edge_fraction"
                    ],
                }
            )
        return result

    selected = last_stage["selected"]
    result.update(
        {
            "selected_stage": last_stage["stage"],
            "matched_guide_count": selected["matched_count"],
            "track_agreement_fraction": selected[
                "agreement_within_4_channels_fraction"
            ],
            "track_correlation": selected["shift_correlation"],
            "left_dynamic_path_z": selected["left_dynamic_path_z"],
            "right_dynamic_path_z": selected["right_dynamic_path_z"],
            "left_search_edge_fraction": selected[
                "left_search_edge_fraction"
            ],
            "right_search_edge_fraction": selected[
                "right_search_edge_fraction"
            ],
        }
    )
    mask_result = evaluate_mask(
        stage_result=last_stage,
        config=config,
        scan=scan,
        beam=beam,
        polar=polar,
        tables_dir=arguments.tables_dir,
        figures_dir=arguments.figures_dir,
        reports_dir=arguments.reports_dir,
        source_dir=source_dir,
        log_path=log_path,
    )
    result["status"] = mask_result["status"]
    result["reason"] = mask_result["reason"]
    result["teeth_above_minimum_snr"] = mask_result.get(
        "teeth_above_minimum_snr", ""
    )
    result["amplitude_snr"] = ";".join(
        f"{value:.6g}" for value in mask_result.get("amplitude_snr", [])
    )
    adaptive = mask_result.get("adaptive_selection")
    if adaptive is not None:
        for key in (
            "half_width_fwhm",
            "masked_tooth_count",
            "masked_fraction",
            "estimated_leakage_fraction",
            "maximum_tooth_leakage_fraction",
        ):
            result[key] = adaptive[key]
    return result


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    summaries = []
    for name, selected in [
        ("all", rows),
        ("XX", [row for row in rows if row["polar"] == "XX"]),
        ("YY", [row for row in rows if row["polar"] == "YY"]),
    ]:
        if not selected:
            continue
        statuses = Counter(str(row["status"]) for row in selected)
        stage_counts = Counter(
            str(row["selected_stage"])
            for row in selected
            if row["selected_stage"]
        )
        usable = statuses["masked"] + statuses["weak_skip"]
        summaries.append(
            {
                "group": name,
                "sample_count": len(selected),
                "masked_count": statuses["masked"],
                "weak_skip_count": statuses["weak_skip"],
                "localization_failed_count": statuses[
                    "localization_failed"
                ],
                "mask_failed_count": statuses["mask_failed"],
                "pipeline_error_count": statuses["pipeline_error"],
                "masked_or_weak_skip_count": usable,
                "masked_or_weak_skip_fraction": usable / len(selected),
                "primary_stage_count": stage_counts["primary"],
                "expanded_stage_count": stage_counts["expanded"],
                "edge_stage_count": stage_counts["edge"],
                "edge_triggered_count": sum(
                    boolean(row["edge_triggered"]) for row in selected
                ),
            }
        )
    return summaries


def main() -> int:
    arguments = parse_arguments()
    with arguments.config.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    selection_rows = read_csv(arguments.selection_csv)
    source_dir = Path(__file__).resolve().parent
    for directory in (
        arguments.tables_dir,
        arguments.figures_dir,
        arguments.reports_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    details = []
    for index, sample in enumerate(arguments.sample, start=1):
        scan, beam, polar = sample
        print(
            f"[{index}/{len(arguments.sample)}] "
            f"scan={scan}, beam={beam}, polar={polar}",
            flush=True,
        )
        try:
            row = process_sample(
                sample=sample,
                selection_rows=selection_rows,
                config=config,
                config_path=arguments.config.resolve(),
                arguments=arguments,
                source_dir=source_dir,
            )
        except Exception as error:
            try:
                date = find_selection(
                    selection_rows, scan, beam
                )["date"]
            except Exception:
                date = ""
            row = empty_result(scan, beam, polar, date)
            row["status"] = "pipeline_error"
            row["reason"] = f"{type(error).__name__}: {error}"
            traceback.print_exc()
            if not arguments.keep_going:
                raise
        details.append(row)
        print(
            f"result={row['status']}, stage={row['selected_stage']}, "
            f"reason={row['reason']}",
            flush=True,
        )
        write_csv(arguments.output_detail_csv, details)

    summaries = summarize(details)
    write_csv(arguments.output_summary_csv, summaries)
    print(f"saved_detail_csv={arguments.output_detail_csv}")
    print(f"saved_summary_csv={arguments.output_summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
