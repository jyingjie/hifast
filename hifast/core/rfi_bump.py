"""Detection and masking functions for quasi-periodic narrow bump RFI."""

from __future__ import annotations

__all__ = [
    "ALGORITHM",
    "CONFIG_PATH",
    "apply_masks",
    "detect_masks",
    "nearest_frequency_indices",
]

import csv
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from .. import rfi_bump_support


ALGORITHM = "hierarchical_track_shape_gate_v5_3"
MINIMUM_AMPLITUDE_SNR = 0.5
MAXIMUM_FWHM_KHZ = 120.0
MARGIN_CHANNELS = 0.5
POLARIZATION_INDEX = {"XX": 0, "YY": 1}
SUPPORT_DIR = Path(rfi_bump_support.__file__).resolve().parent
CONFIG_PATH = (
    SUPPORT_DIR / "configs" / "dynamic_waterfall_mask_v5_3_shape_holdout.yaml"
)


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def prepare_selection(
    fpath,
    destination,
    beam_number,
    date,
    start,
    stop,
):
    write_csv(
        destination,
        [
            {
                "scan": 0,
                "beam": beam_number,
                "date": date,
                "stage_path": str(Path(fpath).resolve()),
                "time_index_start": start,
                "time_index_stop_exclusive": stop,
            }
        ],
    )


def run_locator(
    *,
    fpath,
    beam_number,
    beam,
    date,
    polarizations,
    start,
    stop,
    workdir,
):
    workdir = Path(workdir)
    tables = workdir / "tables"
    figures = workdir / "figures"
    reports = workdir / "reports"
    for directory in (tables, figures, reports):
        directory.mkdir(parents=True, exist_ok=True)
    selection = workdir / "selection.csv"
    prepare_selection(
        fpath,
        selection,
        beam_number,
        date,
        start,
        stop,
    )
    detail = workdir / "detail.csv"
    summary = workdir / "summary.csv"
    command = [
        sys.executable,
        str(SUPPORT_DIR / "run_hierarchical_mask_validation.py"),
        "--config",
        str(CONFIG_PATH),
        "--selection-csv",
        str(selection),
    ]
    for polar in polarizations:
        command.extend(["--sample", f"0|{beam}|{polar}"])
    command.extend(
        [
            "--tables-dir",
            str(tables),
            "--figures-dir",
            str(figures),
            "--reports-dir",
            str(reports),
            "--output-label",
            "rfi_bump_v5_3",
            "--output-detail-csv",
            str(detail),
            "--output-summary-csv",
            str(summary),
            "--no-reuse-waterfalls",
        ]
    )
    environment = os.environ.copy()
    environment.setdefault(
        "MPLCONFIGDIR",
        str(workdir / "matplotlib_cache"),
    )
    completed = subprocess.run(command, env=environment)
    if completed.returncode:
        raise RuntimeError(
            f"bump locator exited with code {completed.returncode}"
        )
    rows = read_csv(detail)
    if len(rows) != len(polarizations):
        raise RuntimeError(
            f"locator returned {len(rows)} rows for "
            f"{len(polarizations)} polarizations"
        )
    pipeline_errors = [
        row for row in rows if row["status"] == "pipeline_error"
    ]
    if pipeline_errors:
        reasons = "; ".join(
            f"{row['polar']}: {row['reason']}"
            for row in pipeline_errors
        )
        raise RuntimeError(
            f"bump locator pipeline error: {reasons}"
        )
    return rows, tables


def group_mask(
    frequency_mhz,
    raw_shift_mhz,
    centers_mhz,
    half_widths_mhz,
):
    mask = np.zeros(
        (len(raw_shift_mhz), len(frequency_mhz)),
        dtype=bool,
    )
    for center_mhz, half_width_mhz in zip(
        centers_mhz,
        half_widths_mhz,
    ):
        dynamic_center = center_mhz + raw_shift_mhz
        mask |= (
            np.abs(
                frequency_mhz[None, :]
                - dynamic_center[:, None]
            )
            <= half_width_mhz
        )
    return mask


def reconstruct_mask(detail, tables, beam):
    stage = detail["selected_stage"]
    polar = detail["polar"]
    candidate = Path(tables) / (
        f"dynamic_mask_scan0_{beam}_{polar}_rfi_bump_v5_3_"
        f"{beam.lower()}_{stage}_mask_eval_candidate_masks.npz"
    )
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    with np.load(candidate, allow_pickle=False) as product:
        frequency_mhz = product["frequency_mhz"].astype(
            np.float64
        )
        original_indices = product[
            "selected_original_indices"
        ].astype(np.int64)
        shift_mhz = product[
            "raw_common_shift_mhz"
        ].astype(np.float64)
        centers_mhz = product[
            "fitted_centers_mhz"
        ].astype(np.float64)
        fwhm_mhz = product[
            "fitted_fwhm_mhz"
        ].astype(np.float64)
    amplitude_snr = np.asarray(
        [
            float(value)
            for value in detail["amplitude_snr"].split(";")
        ],
        dtype=np.float64,
    )
    if len(amplitude_snr) != len(fwhm_mhz):
        raise ValueError(
            f"S/N and fitted-tooth counts differ for {beam}/{polar}"
        )
    selected = (
        (amplitude_snr >= MINIMUM_AMPLITUDE_SNR)
        & (fwhm_mhz * 1000.0 <= MAXIMUM_FWHM_KHZ)
    )
    expected = int(detail["masked_tooth_count"])
    if np.count_nonzero(selected) != expected:
        raise ValueError(
            f"reconstructed tooth count differs for {beam}/{polar}: "
            f"{np.count_nonzero(selected)} != {expected}"
        )
    channel_width_mhz = abs(
        float(np.median(np.diff(frequency_mhz)))
    )
    half_widths_mhz = (
        float(detail["half_width_fwhm"]) * fwhm_mhz[selected]
        + MARGIN_CHANNELS * channel_width_mhz
    )
    mask = group_mask(
        frequency_mhz,
        shift_mhz,
        centers_mhz[selected],
        half_widths_mhz,
    )
    measured = float(np.mean(mask))
    recorded = float(detail["masked_fraction"])
    if not np.isclose(measured, recorded, atol=1e-12):
        raise ValueError(
            f"mask fraction differs for {beam}/{polar}: "
            f"{measured} != {recorded}"
        )
    return original_indices, frequency_mhz, mask


def nearest_frequency_indices(
    source_frequency_mhz,
    target_frequency_mhz,
):
    indices = np.searchsorted(
        source_frequency_mhz,
        target_frequency_mhz,
    )
    indices = np.clip(
        indices,
        0,
        len(source_frequency_mhz) - 1,
    )
    left = np.clip(
        indices - 1,
        0,
        len(source_frequency_mhz) - 1,
    )
    use_left = (
        np.abs(
            source_frequency_mhz[left] - target_frequency_mhz
        )
        < np.abs(
            source_frequency_mhz[indices]
            - target_frequency_mhz
        )
    )
    indices[use_left] = left[use_left]
    difference = np.abs(
        source_frequency_mhz[indices] - target_frequency_mhz
    )
    if float(np.max(difference)) > 1e-6:
        raise ValueError(
            "locator and input frequency grids differ by more than 1 Hz"
        )
    if not np.array_equal(
        indices,
        np.arange(indices[0], indices[0] + len(indices)),
    ):
        raise ValueError(
            "target frequency channels are not contiguous"
        )
    return indices


def build_masks(rows, tables, beam, source_frequency_mhz):
    masks = {}
    results = []
    for row in rows:
        result = {
            "polar": row["polar"],
            "status": row["status"],
            "selected_stage": row["selected_stage"],
            "reason": row["reason"],
            "masked_tooth_count": (
                int(row["masked_tooth_count"])
                if row["masked_tooth_count"]
                else 0
            ),
            "masked_fraction_target_band": (
                float(row["masked_fraction"])
                if row["masked_fraction"]
                else 0.0
            ),
            "masked_sample_count": 0,
            "new_nan_count": 0,
        }
        if row["status"] == "masked":
            time_indices, frequency_mhz, mask = reconstruct_mask(
                row,
                tables,
                beam,
            )
            frequency_indices = nearest_frequency_indices(
                source_frequency_mhz,
                frequency_mhz,
            )
            masks[POLARIZATION_INDEX[row["polar"]]] = (
                time_indices,
                frequency_indices,
                mask,
            )
            result["masked_sample_count"] = int(
                np.count_nonzero(mask)
            )
        results.append(result)
    return masks, results


def detect_masks(
    *,
    fpath,
    beam_number,
    date,
    polarizations,
    start,
    stop,
    workdir,
    source_frequency_mhz,
):
    beam = f"M{beam_number:02d}"
    rows, tables = run_locator(
        fpath=fpath,
        beam_number=beam_number,
        beam=beam,
        date=date,
        polarizations=polarizations,
        start=start,
        stop=stop,
        workdir=workdir,
    )
    return build_masks(
        rows,
        tables,
        beam,
        np.asarray(source_frequency_mhz, dtype=np.float64),
    )


def apply_masks(spectra, masks, results):
    output = np.array(spectra, copy=True)
    result_by_polar = {
        POLARIZATION_INDEX[result["polar"]]: result
        for result in results
    }
    for polar_index, (
        time_indices,
        frequency_indices,
        mask,
    ) in masks.items():
        if not np.array_equal(
            time_indices,
            np.arange(
                time_indices[0],
                time_indices[0] + len(time_indices),
            ),
        ):
            raise ValueError(
                "target time samples are not contiguous"
            )
        time_slice = slice(
            int(time_indices[0]),
            int(time_indices[-1]) + 1,
        )
        frequency_slice = slice(
            int(frequency_indices[0]),
            int(frequency_indices[-1]) + 1,
        )
        block = output[
            time_slice,
            frequency_slice,
            polar_index,
        ]
        if block.shape != mask.shape:
            raise ValueError(
                f"mask shape {mask.shape} differs from "
                f"spectra selection {block.shape}"
            )
        result_by_polar[polar_index]["new_nan_count"] = int(
            np.count_nonzero(mask & np.isfinite(block))
        )
        block[mask] = np.nan
        output[
            time_slice,
            frequency_slice,
            polar_index,
        ] = block
    return output
