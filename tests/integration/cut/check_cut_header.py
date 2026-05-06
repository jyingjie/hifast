#!/usr/bin/env python3
"""
Compare raw FAST FITS headers against a `hifast.cut` HDF5 output.

This script focuses on three checks:
1. Whether `/0.attrs` matches the input FITS primary header.
2. Whether `/1.attrs` mostly matches the input FITS table header, except for
   fields that `hifast.cut` intentionally rewrites after merge/crop.
3. Whether `/Header.attrs` only carries HiFAST history metadata rather than
   duplicating FITS header content.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from collections import defaultdict
from glob import glob
from typing import Dict, Iterable, List, Tuple

import h5py
from astropy.io import fits


EXPECTED_H1_REWRITES = {"NAXIS2", "TDIM21"}


def normalize_value(value):
    """Convert values to JSON-friendly scalars/strings for stable comparison."""
    try:
        import numpy as np  # local import to keep the script standalone
    except ImportError:  # pragma: no cover
        np = None

    if np is not None and isinstance(value, np.generic):
        return value.item()
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except UnicodeDecodeError:
            return repr(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def load_fits_headers(path: str) -> Tuple[Dict[str, object], Dict[str, object]]:
    """Read HDU0 and HDU1 headers without loading the large table data."""
    h0 = dict(fits.getheader(path, 0).items())
    h1 = dict(fits.getheader(path, 1).items())
    return h0, h1


def infer_chunk_paths(first_chunk: str, start: int, stop: int | None) -> List[str]:
    """Infer the FITS chunk list from a `_0001.fits`-style path."""
    base = re.sub(r"[0-9]{4}\.fits\Z", "", first_chunk)
    if base == first_chunk:
        raise ValueError(f"Cannot infer chunk pattern from {first_chunk}")
    if stop is None:
        stop = len(glob(base + "*.fits"))
    return [f"{base}{idx:04d}.fits" for idx in range(start, stop + 1)]


def find_varying_keys(headers: Iterable[Dict[str, object]]) -> Dict[str, List[object]]:
    """Return keys whose values vary across chunks."""
    values = defaultdict(list)
    for header in headers:
        for key, value in header.items():
            values[key].append(normalize_value(value))
    varying = {}
    for key, seq in values.items():
        uniq = []
        for value in seq:
            if value not in uniq:
                uniq.append(value)
        if len(uniq) > 1:
            varying[key] = uniq
    return varying


def compare_headers(expected: Dict[str, object], actual: Dict[str, object]) -> Dict[str, object]:
    """Compare two attribute dictionaries."""
    expected_n = {k: normalize_value(v) for k, v in expected.items()}
    actual_n = {k: normalize_value(v) for k, v in actual.items()}

    missing = sorted(k for k in expected_n.keys() if k not in actual_n)
    extra = sorted(k for k in actual_n.keys() if k not in expected_n)
    differing = {}
    for key in sorted(expected_n.keys() & actual_n.keys()):
        if expected_n[key] != actual_n[key]:
            differing[key] = {"expected": expected_n[key], "actual": actual_n[key]}
    return {"missing_keys": missing, "extra_keys": extra, "differing_values": differing}


def load_hdf5_summary(path: str) -> Dict[str, object]:
    """Read the HDF5 groups/attrs needed for header comparison."""
    with h5py.File(path, "r") as f:
        g0 = dict(f["0"].attrs.items()) if "0" in f else {}
        g1 = dict(f["1"].attrs.items()) if "1" in f else {}
        gh = dict(f["Header"].attrs.items()) if "Header" in f else {}

        data_shape = tuple(int(v) for v in f["1"]["DATA"].shape) if "1" in f and "DATA" in f["1"] else None
        nchan_col = int(f["1"]["NCHAN"][0]) if "1" in f and "NCHAN" in f["1"] else None
        freq_col = normalize_value(f["1"]["FREQ"][0]) if "1" in f and "FREQ" in f["1"] else None
        row_count = int(f["1"]["UTOBS"].shape[0]) if "1" in f and "UTOBS" in f["1"] else None

    return {
        "group0_attrs": g0,
        "group1_attrs": g1,
        "header_attrs": gh,
        "data_shape": data_shape,
        "nchan_column_first": nchan_col,
        "freq_column_first": freq_col,
        "row_count": row_count,
    }


def build_expected_rewrites(hdf5_info: Dict[str, object]) -> Dict[str, object]:
    """Fields that `hifast.cut` intentionally updates in `/1.attrs`."""
    row_count = hdf5_info["row_count"]
    nchan = hdf5_info["nchan_column_first"]
    return {
        "NAXIS2": row_count,
        "TDIM21": f"(2, {nchan})" if nchan is not None else None,
    }


def split_h1_differences(diff: Dict[str, object], expected_rewrites: Dict[str, object]) -> Tuple[Dict[str, object], Dict[str, object]]:
    """
    Split `/1` differences into expected rewrites and unexpected differences.
    """
    differing = diff["differing_values"]
    expected = {}
    unexpected = {}
    for key, payload in differing.items():
        expected_value = normalize_value(expected_rewrites.get(key))
        if key in EXPECTED_H1_REWRITES and expected_value == normalize_value(payload["actual"]):
            expected[key] = payload
        else:
            unexpected[key] = payload
    return expected, unexpected


def summarize_header_overlap(fits_h0: Dict[str, object], fits_h1: Dict[str, object], history_header: Dict[str, object]) -> Dict[str, List[str]]:
    """
    Report overlap between HiFAST history attrs and FITS header keys.
    """
    history_keys = set(history_header.keys())
    return {
        "overlap_with_hdu0": sorted(history_keys & set(fits_h0.keys())),
        "overlap_with_hdu1": sorted(history_keys & set(fits_h1.keys())),
    }


def build_report(first_chunk: str, hdf5_path: str, start: int, stop: int | None) -> Dict[str, object]:
    chunk_paths = infer_chunk_paths(first_chunk, start, stop)
    h0_list = []
    h1_list = []
    for path in chunk_paths:
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        h0, h1 = load_fits_headers(path)
        h0_list.append(h0)
        h1_list.append(h1)

    ref_h0 = h0_list[0]
    ref_h1 = h1_list[0]
    hdf5_info = load_hdf5_summary(hdf5_path)
    h0_diff = compare_headers(ref_h0, hdf5_info["group0_attrs"])
    h1_diff = compare_headers(ref_h1, hdf5_info["group1_attrs"])
    expected_rewrites = build_expected_rewrites(hdf5_info)
    expected_h1_diff, unexpected_h1_diff = split_h1_differences(h1_diff, expected_rewrites)
    header_overlap = summarize_header_overlap(ref_h0, ref_h1, hdf5_info["header_attrs"])

    return {
        "inputs": {
            "first_chunk": first_chunk,
            "chunk_count": len(chunk_paths),
            "chunk_paths": chunk_paths,
            "hdf5_path": hdf5_path,
        },
        "chunk_header_variation": {
            "hdu0_varying_keys": find_varying_keys(h0_list),
            "hdu1_varying_keys": find_varying_keys(h1_list),
        },
        "hdf5_summary": {
            "data_shape": hdf5_info["data_shape"],
            "row_count": hdf5_info["row_count"],
            "nchan_column_first": hdf5_info["nchan_column_first"],
            "freq_column_first": hdf5_info["freq_column_first"],
        },
        "group0_vs_fits_hdu0": h0_diff,
        "group1_vs_fits_hdu1": {
            "missing_keys": h1_diff["missing_keys"],
            "extra_keys": h1_diff["extra_keys"],
            "expected_rewritten_values": expected_h1_diff,
            "unexpected_differing_values": unexpected_h1_diff,
        },
        "expected_group1_rewrites": expected_rewrites,
        "history_header_overlap": header_overlap,
        "history_header_keys": sorted(hdf5_info["header_attrs"].keys()),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check whether hifast.cut preserved FITS headers correctly."
    )
    parser.add_argument("first_chunk", help="Path to the first FITS chunk, e.g. *_0001.fits")
    parser.add_argument("hdf5_path", help="Path to the HDF5 file produced by hifast.cut")
    parser.add_argument("--start", type=int, default=1, help="First chunk index used by hifast.cut")
    parser.add_argument("--stop", type=int, help="Last chunk index used by hifast.cut")
    parser.add_argument(
        "--json-out",
        help="Optional path to save the full JSON report",
    )
    parser.add_argument(
        "--fail-on-unexpected",
        action="store_true",
        help="Exit with code 1 if unexpected header differences are found.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build_report(args.first_chunk, args.hdf5_path, args.start, args.stop)

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, ensure_ascii=False)

    print(json.dumps(report, indent=2, ensure_ascii=False))
    if args.fail_on_unexpected:
        if report["group0_vs_fits_hdu0"]["missing_keys"]:
            return 1
        if report["group0_vs_fits_hdu0"]["extra_keys"]:
            return 1
        if report["group0_vs_fits_hdu0"]["differing_values"]:
            return 1
        if report["group1_vs_fits_hdu1"]["missing_keys"]:
            return 1
        if report["group1_vs_fits_hdu1"]["extra_keys"]:
            return 1
        if report["group1_vs_fits_hdu1"]["unexpected_differing_values"]:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
