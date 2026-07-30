#!/usr/bin/env python3
"""Fit de-drifted Gaussian profiles for both interleaved comb groups."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import shift as shift_row
from scipy.optimize import curve_fit


def gaussian_with_linear_background(
    frequency_mhz: np.ndarray,
    constant: float,
    slope: float,
    amplitude: float,
    center_mhz: float,
    sigma_mhz: float,
) -> np.ndarray:
    return (
        constant
        + slope * (frequency_mhz - center_mhz)
        + amplitude
        * np.exp(
            -0.5 * ((frequency_mhz - center_mhz) / sigma_mhz) ** 2
        )
    )


def robust_scale(values: np.ndarray) -> float:
    center = float(np.nanmedian(values))
    scale = float(1.4826 * np.nanmedian(np.abs(values - center)))
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanstd(values))
    return scale


def fit_profile(
    frequency_mhz: np.ndarray,
    profile: np.ndarray,
    guide_center_mhz: float,
    fit_half_width_mhz: float,
    center_bound_mhz: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    selected = np.abs(frequency_mhz - guide_center_mhz) <= fit_half_width_mhz
    x = frequency_mhz[selected]
    y = profile[selected]
    finite = np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if len(x) < 12:
        raise RuntimeError("Too few finite channels for a profile fit")
    central = np.abs(x - guide_center_mhz) <= center_bound_mhz
    baseline = float(np.median(y))
    peak_index = int(np.argmax(y[central]))
    initial_center = float(x[central][peak_index])
    initial_amplitude = max(
        float(np.max(y[central]) - baseline), np.finfo(float).eps
    )
    fitted, covariance = curve_fit(
        gaussian_with_linear_background,
        x,
        y,
        p0=[baseline, 0.0, initial_amplitude, initial_center, 0.025],
        bounds=(
            [
                -np.inf,
                -np.inf,
                0.0,
                guide_center_mhz - center_bound_mhz,
                0.006,
            ],
            [
                np.inf,
                np.inf,
                np.inf,
                guide_center_mhz + center_bound_mhz,
                0.08,
            ],
        ),
        maxfev=5000,
    )
    residual = y - gaussian_with_linear_background(x, *fitted)
    noise = robust_scale(residual)
    return fitted, covariance, noise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--waterfall", type=Path, required=True)
    parser.add_argument("--group-a", type=Path, required=True)
    parser.add_argument("--group-b", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-npz", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--fit-half-width-mhz", type=float, default=0.12)
    parser.add_argument("--center-bound-mhz", type=float, default=0.05)
    args = parser.parse_args()

    with np.load(args.waterfall, allow_pickle=False) as product:
        frequency_mhz = product["frequency_mhz"].astype(np.float64)
        binned_residual = product["binned_residual"].astype(np.float64)
        data_unit = str(product["data_unit"])
    with np.load(args.group_a, allow_pickle=False) as product:
        guides_a = product["guide_centers_mhz"].astype(np.float64)
        common_shift_channels = product[
            "common_shift_channels"
        ].astype(np.float64)
    with np.load(args.group_b, allow_pickle=False) as product:
        guides_b = product["guide_centers_mhz"].astype(np.float64)

    aligned = np.empty_like(binned_residual)
    for time_index, (row, shift_channels) in enumerate(
        zip(binned_residual, common_shift_channels)
    ):
        aligned[time_index] = shift_row(
            row,
            -shift_channels,
            order=1,
            mode="nearest",
            prefilter=False,
        )
    stacked_profile = np.nanmean(aligned, axis=0)

    rows: list[dict[str, object]] = []
    fitted_centers = []
    fitted_fwhm = []
    fitted_amplitudes = []
    fitted_parameters = []
    fitted_groups = []
    fitted_tooth_numbers = []
    fitted_models: list[tuple[np.ndarray, np.ndarray, str]] = []
    for group, guides in (("A_stronger", guides_a), ("B_weaker", guides_b)):
        for tooth_number, guide_center in enumerate(guides, start=1):
            fitted, covariance, noise = fit_profile(
                frequency_mhz,
                stacked_profile,
                float(guide_center),
                args.fit_half_width_mhz,
                args.center_bound_mhz,
            )
            center_mhz = float(fitted[3])
            sigma_mhz = float(fitted[4])
            fwhm_mhz = 2.354820045 * sigma_mhz
            amplitude = float(fitted[2])
            center_uncertainty_mhz = float(
                np.sqrt(max(0.0, covariance[3, 3]))
            )
            rows.append(
                {
                    "group": group,
                    "tooth_number": tooth_number,
                    "guide_center_mhz": f"{guide_center:.9f}",
                    "fitted_center_mhz": f"{center_mhz:.9f}",
                    "center_offset_khz": (
                        f"{(center_mhz - guide_center) * 1000:.6f}"
                    ),
                    "center_uncertainty_khz": (
                        f"{center_uncertainty_mhz * 1000:.6f}"
                    ),
                    "fwhm_khz": f"{fwhm_mhz * 1000:.6f}",
                    "amplitude": f"{amplitude:.12g}",
                    "residual_noise": f"{noise:.12g}",
                    "amplitude_snr": (
                        f"{amplitude / noise:.6f}"
                        if noise > 0
                        else "nan"
                    ),
                    "data_unit": data_unit,
                }
            )
            fitted_centers.append(center_mhz)
            fitted_fwhm.append(fwhm_mhz)
            fitted_amplitudes.append(amplitude)
            fitted_parameters.append(fitted)
            fitted_groups.append(group)
            fitted_tooth_numbers.append(tooth_number)
            selected = (
                np.abs(frequency_mhz - guide_center)
                <= args.fit_half_width_mhz
            )
            fitted_models.append(
                (
                    frequency_mhz[selected],
                    gaussian_with_linear_background(
                        frequency_mhz[selected], *fitted
                    ),
                    group,
                )
            )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        frequency_mhz=frequency_mhz,
        stacked_profile=stacked_profile,
        fitted_groups=np.asarray(fitted_groups),
        fitted_tooth_numbers=np.asarray(fitted_tooth_numbers, dtype=np.int64),
        fitted_centers_mhz=np.asarray(fitted_centers),
        fitted_fwhm_mhz=np.asarray(fitted_fwhm),
        fitted_amplitudes=np.asarray(fitted_amplitudes),
        fitted_parameters=np.asarray(fitted_parameters),
        common_shift_source=str(args.group_a),
        source_waterfall=str(args.waterfall),
        data_unit=data_unit,
    )

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(15, 8),
        gridspec_kw={"height_ratios": [1.8, 0.8]},
    )
    scale = 1000 if data_unit in {"Jy/beam", "K"} else 1
    display_unit = (
        "mJy/beam"
        if data_unit == "Jy/beam"
        else "mK"
        if data_unit == "K"
        else data_unit
    )
    axes[0].plot(
        frequency_mhz,
        stacked_profile * scale,
        color="0.35",
        lw=0.8,
        label="de-drifted mean residual",
    )
    colors = {"A_stronger": "black", "B_weaker": "tab:green"}
    labels_used: set[str] = set()
    for x, model, group in fitted_models:
        label = group if group not in labels_used else None
        labels_used.add(group)
        axes[0].plot(
            x,
            model * scale,
            color=colors[group],
            lw=1.3,
            label=label,
        )
    axes[0].set_ylabel(f"residual ({display_unit})")
    axes[0].set_xlabel("de-drifted topocentric frequency (MHz)")
    axes[0].set_title("Independent Gaussian profile fits after A-shift removal")
    axes[0].grid(alpha=0.2)
    axes[0].legend()

    x_positions = np.arange(len(rows))
    bar_colors = [colors[str(row["group"])] for row in rows]
    axes[1].bar(
        x_positions,
        np.asarray(fitted_fwhm) * 1000,
        color=bar_colors,
        alpha=0.75,
    )
    axes[1].set_xticks(
        x_positions,
        [
            f"{row['group'][0]}{row['tooth_number']}"
            for row in rows
        ],
    )
    axes[1].set_ylabel("fitted FWHM (kHz)")
    axes[1].set_xlabel("comb group and tooth")
    axes[1].grid(axis="y", alpha=0.2)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=170)
    plt.close(figure)

    print(f"fitted_profiles={len(rows)}")
    for group in ("A_stronger", "B_weaker"):
        widths = [
            float(row["fwhm_khz"])
            for row in rows
            if row["group"] == group
        ]
        print(f"{group}_median_fwhm_khz={np.median(widths):.6f}")
    print(f"saved_csv={args.output_csv}")
    print(f"saved_npz={args.output_npz}")
    print(f"saved_figure={args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
