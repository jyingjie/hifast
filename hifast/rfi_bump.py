"""Experimental mask for quasi-periodic narrow bump RFI."""

from __future__ import annotations

__all__ = ["IO", "parser"]

from .utils.io import *

import csv
import json
import subprocess
import tempfile
from copy import deepcopy
from pathlib import Path

import numpy as np


ALGORITHM = "hierarchical_track_shape_gate_v5_3"
FREQUENCY_RANGE_MHZ = (1421.0, 1424.8)
MINIMUM_AMPLITUDE_SNR = 0.5
MAXIMUM_FWHM_KHZ = 120.0
MARGIN_CHANNELS = 0.5
POLARIZATION_INDEX = {"XX": 0, "YY": 1}
SUPPORT_DIR = Path(__file__).with_name("rfi_bump_support")
CONFIG_PATH = (
    SUPPORT_DIR / "configs" / "dynamic_waterfall_mask_v5_3_shape_holdout.yaml"
)


parser = ArgumentParser(
    prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}",
    formatter_class=formatter_class,
    allow_abbrev=False,
    description=(
        "Experimental pre-Doppler mask for the approximately 1 MHz "
        "quasi-periodic narrow bump RFI."
    ),
)
add_common_argument(parser)
add_h5_compression_arguments(parser)
parser.add_argument(
    "fpath",
    help="input pre-Doppler spectra file containing S/flux",
)
parser.add_argument(
    "--polar",
    choices=["both", "XX", "YY"],
    default="both",
    help="polarization to process",
)
parser.add_argument(
    "--time_index_start",
    type=int,
    default=0,
    help="first input spectrum to process",
)
parser.add_argument(
    "--time_index_stop",
    type=int,
    help="exclusive stop index; default: all spectra",
)
parser.add_argument(
    "--keep_products",
    type=bool_fun,
    choices=[True, False],
    default="False",
    help="keep intermediate tables, figures, and logs",
)
parser.add_argument(
    "--products_dir",
    help=(
        "intermediate-product directory; setting it also enables "
        "keep_products"
    ),
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
        "MPLCONFIGDIR", str(workdir / "matplotlib_cache")
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
        raise RuntimeError(f"bump locator pipeline error: {reasons}")
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
        centers_mhz, half_widths_mhz
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
            source_frequency_mhz[indices] - target_frequency_mhz
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


class IO(BaseIO):
    ver = "old"

    def __init__(
        self,
        args,
        dict_in=None,
        inplace_args=False,
        HistoryAdd=None,
    ):
        args = args if inplace_args else deepcopy(args)
        args.no_radec = True
        super().__init__(
            args,
            dict_in=dict_in,
            inplace_args=True,
            HistoryAdd=HistoryAdd,
        )
        if self.infield != "flux":
            raise ValueError(
                "rfi_bump requires an input file containing S/flux"
            )

    def _get_fpart(self):
        return "-rfi_bump"

    def _time_range(self):
        start = self.args.time_index_start
        stop = (
            len(self.mjd)
            if self.args.time_index_stop is None
            else self.args.time_index_stop
        )
        if start < 0 or stop > len(self.mjd) or stop <= start:
            raise ValueError(
                f"invalid time range [{start}, {stop}) for "
                f"{len(self.mjd)} spectra"
            )
        return start, stop

    def _products_path(self):
        args = self.args
        if args.products_dir is not None:
            path = sub_patten(
                args.products_dir,
                date=self.date,
                nB=f"{self.nB:02d}",
                project=self.project,
            )
            return Path(os.path.expanduser(path)).resolve()
        return Path(self.fpath_out).with_suffix("").with_name(
            f"{Path(self.fpath_out).stem}_products"
        )

    def _locate_masks(self, workdir):
        start, stop = self._time_range()
        polarizations = (
            ["XX", "YY"]
            if self.args.polar == "both"
            else [self.args.polar]
        )
        rows, tables = run_locator(
            fpath=self.args.fpath,
            beam_number=self.nB,
            beam=f"M{self.nB:02d}",
            date=self.date,
            polarizations=polarizations,
            start=start,
            stop=stop,
            workdir=workdir,
        )
        return build_masks(
            rows,
            tables,
            f"M{self.nB:02d}",
            np.asarray(self.freq, dtype=np.float64),
        )

    def gen_s2p_out(self):
        keep_products = (
            self.args.keep_products
            or self.args.products_dir is not None
        )
        if keep_products:
            products = self._products_path()
            products.mkdir(parents=True, exist_ok=True)
            masks, results = self._locate_masks(products)
        else:
            with tempfile.TemporaryDirectory(
                prefix=".rfi_bump_",
                dir=self.args.outdir,
            ) as temporary:
                masks, results = self._locate_masks(temporary)
        self.s2p_out = apply_masks(
            self.s2p,
            masks,
            results,
        )
        self.rfi_bump_results = results
        self.Header["rfi_bump_algorithm"] = ALGORITHM
        self.Header["rfi_bump_experimental"] = True
        self.Header["rfi_bump_formal_science_mask"] = False
        self.Header["rfi_bump_results"] = json.dumps(
            results,
            ensure_ascii=False,
        )
        for result in results:
            print(
                f"{result['polar']}: "
                f"status={result['status']}, "
                "masked_samples="
                f"{result['masked_sample_count']}, "
                f"new_nan={result['new_nan_count']}"
            )


if __name__ == "__main__":
    args_ = parser.parse_args()
    print("#" * 35 + "Args" + "#" * 35)
    args_from = parser.format_values()
    print(args_from)
    print("#" * 35 + "####" + "#" * 35)
    HistoryAdd = (
        {"args_from": args_from}
        if args_.my_config is not None
        else None
    )
    io = IO(args_, HistoryAdd=HistoryAdd)
    io()
