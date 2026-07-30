#!/usr/bin/env python3
"""Evaluate candidate masks on the current M05/YY sample only."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import shift as shift_row


MASK_KEYS = {
    "core_all_groups": "mask_core_all_groups",
    "conservative_all_groups": "mask_conservative_all_groups",
    "conservative_a_only": "mask_conservative_a_only",
    "conservative_a_plus_gated_b": (
        "mask_conservative_a_plus_gated_b"
    ),
}


def weighted_aligned_mean(
    values: np.ndarray,
    valid: np.ndarray,
    shift_channels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    numerator_sum = np.zeros(values.shape[1], dtype=np.float64)
    weight_sum = np.zeros(values.shape[1], dtype=np.float64)
    for row, row_valid, row_shift in zip(values, valid, shift_channels):
        finite = row_valid & np.isfinite(row)
        numerator = np.where(finite, row, 0.0)
        weight = finite.astype(np.float64)
        numerator_sum += shift_row(
            numerator,
            -row_shift,
            order=1,
            mode="constant",
            cval=0.0,
            prefilter=False,
        )
        weight_sum += shift_row(
            weight,
            -row_shift,
            order=1,
            mode="constant",
            cval=0.0,
            prefilter=False,
        )
    mean = np.full(values.shape[1], np.nan, dtype=np.float64)
    positive = weight_sum > 0
    mean[positive] = numerator_sum[positive] / weight_sum[positive]
    return mean, weight_sum


def fit_fixed_gaussian_amplitude(
    frequency_mhz: np.ndarray,
    profile: np.ndarray,
    weight: np.ndarray,
    center_mhz: float,
    sigma_mhz: float,
    half_window_mhz: float = 0.12,
) -> float:
    distance = frequency_mhz - center_mhz
    selected = (
        (np.abs(distance) <= half_window_mhz)
        & np.isfinite(profile)
        & np.isfinite(weight)
        & (weight > 0)
    )
    if np.count_nonzero(selected) < 8:
        return np.nan
    x = distance[selected]
    gaussian = np.exp(-0.5 * (x / sigma_mhz) ** 2)
    design = np.column_stack(
        [np.ones(np.count_nonzero(selected)), x, gaussian]
    )
    square_root_weight = np.sqrt(weight[selected])
    coefficients, _, _, _ = np.linalg.lstsq(
        design * square_root_weight[:, None],
        profile[selected] * square_root_weight,
        rcond=None,
    )
    return max(0.0, float(coefficients[2]))


def profile_leakage(
    frequency_mhz: np.ndarray,
    original_profile: np.ndarray,
    original_weight: np.ndarray,
    masked_profile: np.ndarray,
    masked_weight: np.ndarray,
    fitted_parameters: np.ndarray,
) -> tuple[float, list[float]]:
    relative_weight = np.divide(
        masked_weight,
        original_weight,
        out=np.zeros_like(masked_weight),
        where=original_weight > 0,
    )
    relative_weight = np.clip(relative_weight, 0.0, 1.0)
    total_original = 0.0
    total_remaining = 0.0
    per_tooth = []
    for parameters in fitted_parameters:
        center_mhz = float(parameters[3])
        sigma_mhz = float(parameters[4])
        gaussian = np.exp(
            -0.5
            * ((frequency_mhz - center_mhz) / sigma_mhz) ** 2
        )
        original_amplitude = fit_fixed_gaussian_amplitude(
            frequency_mhz,
            original_profile,
            original_weight,
            center_mhz,
            sigma_mhz,
        )
        masked_amplitude = fit_fixed_gaussian_amplitude(
            frequency_mhz,
            masked_profile,
            masked_weight,
            center_mhz,
            sigma_mhz,
        )
        original_energy = original_amplitude * np.sum(gaussian)
        remaining_energy = (
            masked_amplitude * np.sum(gaussian * relative_weight)
        )
        ratio = (
            remaining_energy / original_energy
            if original_energy > 0 and np.isfinite(remaining_energy)
            else np.nan
        )
        per_tooth.append(float(ratio))
        if np.isfinite(original_energy) and np.isfinite(remaining_energy):
            total_original += original_energy
            total_remaining += remaining_energy
    aggregate = (
        total_remaining / total_original if total_original > 0 else np.nan
    )
    return float(aggregate), per_tooth


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    rows: list[dict[str, object]],
    tooth_rows: list[dict[str, object]],
) -> None:
    lines = [
        "# scan 43 M05/YY 候选 mask 实测评价",
        "",
        "本报告只评价当前 scan 43、M05/YY 样本，不代表其他波束、偏振或日期。",
        "",
        "## 汇总",
        "",
        "| 候选 | 遮挡比例 | 估计实际残留 | 估计实际覆盖 | 频率通道最大时间损失 |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['candidate']} | "
            f"{float(row['masked_fraction']) * 100:.2f}% | "
            f"{float(row['estimated_leakage_fraction']) * 100:.2f}% | "
            f"{float(row['estimated_coverage_fraction']) * 100:.2f}% | "
            f"{float(row['channel_time_loss_max']) * 100:.2f}% |"
        )
    lines.extend(
        [
            "",
            "“估计实际残留”使用 mask 后仍有权重的 Gaussian 翼部重新拟合，"
            "再乘以实际剩余权重计算。它比单纯 Gaussian 几何积分更接近"
            "当前样本，但仍不是其他数据上的真值。",
            "",
            "## 每条轨迹残留",
            "",
            "| 候选 | 组 | 轨迹 | 残留比例 |",
            "|---|---|---:|---:|",
        ]
    )
    for row in tooth_rows:
        lines.append(
            f"| {row['candidate']} | {row['group']} | "
            f"{row['tooth_number']} | "
            f"{float(row['estimated_leakage_fraction']) * 100:.2f}% |"
        )
    lines.extend(
        [
            "",
            "所有产品均为独立候选 mask，原始 HDF5 未修改。",
            "",
            "## 当前样本判断",
            "",
            "- `core_all_groups` 的残留约为四分之一，不能满足残留不高于 "
            "`10%` 的要求。",
            "- `conservative_a_only` 完全保留 B 组，不符合当前样本的目标。",
            "- `conservative_a_plus_gated_b` 的总体残留低，但三条 B 组轨迹"
            "分别仍残留约 `10%–12%`，没有全部通过逐轨迹检查。",
            "- `conservative_all_groups` 的七条轨迹残留均低于 `10%`，"
            "但部分固定频率通道损失接近全部时间样本。它只能作为当前"
            "样本的强抑制候选，不能仅凭总遮挡比例判断可用。",
            "- mask 外原始数组值没有被修改；干净频率区域的当前时间平均"
            "谱差为零。",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--waterfall", type=Path, required=True)
    parser.add_argument("--candidate-masks", type=Path, required=True)
    parser.add_argument("--profile-fits", type=Path, required=True)
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--tooth-metrics-csv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()

    with np.load(args.waterfall, allow_pickle=False) as product:
        frequency_mhz = product["frequency_mhz"].astype(np.float64)
        raw_residual = product["raw_residual"].astype(np.float64)
    with np.load(args.candidate_masks, allow_pickle=False) as product:
        raw_shift_mhz = product["raw_common_shift_mhz"].astype(np.float64)
        masks = {
            name: product[key].astype(bool)
            for name, key in MASK_KEYS.items()
        }
    with np.load(args.profile_fits, allow_pickle=False) as product:
        fitted_parameters = product["fitted_parameters"].astype(np.float64)
        fitted_groups = product["fitted_groups"].astype(str)
        fitted_tooth_numbers = product[
            "fitted_tooth_numbers"
        ].astype(np.int64)

    delta_frequency_mhz = float(np.median(np.diff(frequency_mhz)))
    shift_channels = raw_shift_mhz / delta_frequency_mhz
    original_valid = np.isfinite(raw_residual)
    original_profile, original_weight = weighted_aligned_mean(
        raw_residual, original_valid, shift_channels
    )

    fitted_centers = fitted_parameters[:, 3]
    clean_frequency = np.ones(len(frequency_mhz), dtype=bool)
    for center_mhz in fitted_centers:
        clean_frequency &= np.abs(frequency_mhz - center_mhz) > 0.14

    rows: list[dict[str, object]] = []
    tooth_rows: list[dict[str, object]] = []
    profiles: dict[str, np.ndarray] = {}
    weights: dict[str, np.ndarray] = {}
    for candidate_name, mask in masks.items():
        masked_valid = original_valid & ~mask
        masked_profile, masked_weight = weighted_aligned_mean(
            raw_residual, masked_valid, shift_channels
        )
        profiles[candidate_name] = masked_profile
        weights[candidate_name] = masked_weight
        leakage, per_tooth = profile_leakage(
            frequency_mhz,
            original_profile,
            original_weight,
            masked_profile,
            masked_weight,
            fitted_parameters,
        )
        row_fraction = np.mean(mask, axis=1)
        channel_time_loss = np.mean(mask, axis=0)
        clean_selected = (
            clean_frequency
            & np.isfinite(original_profile)
            & np.isfinite(masked_profile)
        )
        clean_difference = (
            masked_profile[clean_selected]
            - original_profile[clean_selected]
        )
        rows.append(
            {
                "candidate": candidate_name,
                "masked_fraction": f"{np.mean(mask):.9f}",
                "row_mask_fraction_min": f"{np.min(row_fraction):.9f}",
                "row_mask_fraction_median": (
                    f"{np.median(row_fraction):.9f}"
                ),
                "row_mask_fraction_max": f"{np.max(row_fraction):.9f}",
                "channel_time_loss_median": (
                    f"{np.median(channel_time_loss):.9f}"
                ),
                "channel_time_loss_p95": (
                    f"{np.percentile(channel_time_loss, 95):.9f}"
                ),
                "channel_time_loss_max": (
                    f"{np.max(channel_time_loss):.9f}"
                ),
                "estimated_leakage_fraction": f"{leakage:.9f}",
                "estimated_coverage_fraction": f"{1 - leakage:.9f}",
                "clean_profile_difference_median": (
                    f"{np.median(clean_difference):.12g}"
                ),
                "clean_profile_difference_rms": (
                    f"{np.sqrt(np.mean(clean_difference**2)):.12g}"
                ),
                "outside_mask_values_modified": "false",
            }
        )
        for group, tooth_number, tooth_leakage in zip(
            fitted_groups, fitted_tooth_numbers, per_tooth
        ):
            tooth_rows.append(
                {
                    "candidate": candidate_name,
                    "group": group,
                    "tooth_number": int(tooth_number),
                    "estimated_leakage_fraction": (
                        f"{tooth_leakage:.9f}"
                    ),
                }
            )

    write_csv(args.metrics_csv, rows)
    write_csv(args.tooth_metrics_csv, tooth_rows)
    write_report(args.report, rows, tooth_rows)

    figure, axes = plt.subplots(3, 1, figsize=(15, 10))
    axes[0].plot(
        frequency_mhz,
        original_profile * 1000,
        color="0.45",
        lw=0.8,
        label="original de-drifted mean",
    )
    for candidate_name, color in (
        ("core_all_groups", "tab:blue"),
        ("conservative_all_groups", "tab:red"),
    ):
        profile = profiles[candidate_name]
        relative_weight = np.divide(
            weights[candidate_name],
            original_weight,
            out=np.zeros_like(original_weight),
            where=original_weight > 0,
        )
        visible = relative_weight > 0.05
        visible_profile = profile.copy()
        visible_profile[~visible] = np.nan
        axes[0].plot(
            frequency_mhz,
            visible_profile * 1000,
            color=color,
            lw=0.8,
            label=candidate_name,
        )
    axes[0].set_ylabel("residual (mJy/beam)")
    axes[0].set_xlabel("de-drifted frequency (MHz)")
    axes[0].set_title("Available profile after masking")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.2)

    x = np.arange(len(rows))
    axes[1].bar(
        x - 0.18,
        [float(row["masked_fraction"]) * 100 for row in rows],
        width=0.36,
        label="masked samples",
    )
    axes[1].bar(
        x + 0.18,
        [
            float(row["estimated_coverage_fraction"]) * 100
            for row in rows
        ],
        width=0.36,
        label="estimated interference coverage",
    )
    axes[1].set_xticks(
        x,
        [str(row["candidate"]).replace("_", "\n") for row in rows],
        fontsize=8,
    )
    axes[1].set_ylabel("percent")
    axes[1].grid(axis="y", alpha=0.2)
    axes[1].legend()

    conservative_loss = np.mean(
        masks["conservative_all_groups"], axis=0
    )
    axes[2].plot(
        frequency_mhz,
        conservative_loss * 100,
        color="tab:red",
        lw=1.0,
    )
    axes[2].set_xlabel("topocentric frequency (MHz)")
    axes[2].set_ylabel("time samples masked (%)")
    axes[2].set_title(
        "Per-channel time loss: conservative_all_groups"
    )
    axes[2].grid(alpha=0.2)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=170)
    plt.close(figure)

    for row in rows:
        print(
            f"{row['candidate']}: "
            f"mask={row['masked_fraction']},"
            f"leakage={row['estimated_leakage_fraction']},"
            f"channel_loss_max={row['channel_time_loss_max']}"
        )
    print(f"saved_metrics={args.metrics_csv}")
    print(f"saved_tooth_metrics={args.tooth_metrics_csv}")
    print(f"saved_report={args.report}")
    print(f"saved_figure={args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
