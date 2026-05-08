#!/usr/bin/env python3
"""
Benchmark compression choices for `hifast.cut` outputs on real HDF5 data.

The script rewrites `/1/DATA` from an existing cut output into standalone HDF5
files using different compression/chunk settings, then reports file size and
simple read/write timings.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import h5py

try:
    import hdf5plugin
except ImportError:  # pragma: no cover
    hdf5plugin = None


DEFAULT_CASES = [
    {"name": "gzip_2_chunk_2x64x512", "kind": "gzip", "level": 2, "shuffle": True, "chunk": (2, 64, 512)},
    {"name": "gzip_3_chunk_2x64x1024", "kind": "gzip", "level": 3, "shuffle": True, "chunk": (2, 64, 1024)},
    {"name": "gzip_7_chunk_2x128x512", "kind": "gzip", "level": 7, "shuffle": True, "chunk": (2, 128, 512)},
    {"name": "gzip_9_chunk_2x64x1024", "kind": "gzip", "level": 9, "shuffle": True, "chunk": (2, 64, 1024)},
    {"name": "lzf_chunk_2x64x1024", "kind": "lzf", "shuffle": True, "chunk": (2, 64, 1024)},
    {"name": "blosc2_lz4_5_chunk_2x128x1024", "kind": "blosc2_lz4", "level": 5, "chunk": (2, 128, 1024)},
    {"name": "blosc2_zstd_5_chunk_2x128x1024", "kind": "blosc2_zstd", "level": 5, "chunk": (2, 128, 1024)},
    {"name": "bitshuffle_lz4_chunk_2x64x1024", "kind": "bitshuffle_lz4", "chunk": (2, 64, 1024)},
    {"name": "zstd_9_chunk_2x64x1024", "kind": "zstd", "level": 9, "chunk": (2, 64, 1024)},
    {"name": "bitshuffle_zstd_5_chunk_2x64x1024", "kind": "bitshuffle_zstd", "level": 5, "chunk": (2, 64, 1024)},
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark compression for hifast cut outputs.")
    parser.add_argument("source_hdf5", help="Existing cut output HDF5 file.")
    parser.add_argument("--outdir", required=True, help="Directory to store benchmark outputs.")
    parser.add_argument(
        "--case",
        action="append",
        dest="cases",
        help="Case name to run. Can be specified multiple times. Default runs all cases.",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        help="Optional row limit for faster benchmarking on a subset.",
    )
    parser.add_argument(
        "--max-chans",
        type=int,
        help="Optional channel limit for faster benchmarking on a subset.",
    )
    parser.add_argument(
        "--json-out",
        help="Optional path to save the benchmark summary JSON.",
    )
    return parser.parse_args()


def select_cases(case_names: list[str] | None) -> list[dict]:
    if not case_names:
        return DEFAULT_CASES
    name_set = set(case_names)
    selected = [case for case in DEFAULT_CASES if case["name"] in name_set]
    missing = sorted(name_set - {case["name"] for case in selected})
    if missing:
        raise ValueError(f"Unknown case names: {', '.join(missing)}")
    return selected


def resolve_case_kwargs(case: dict) -> dict:
    kind = case["kind"]
    if kind == "gzip":
        return {
            "compression": "gzip",
            "compression_opts": case["level"],
            "shuffle": case.get("shuffle", False),
        }
    if kind == "lzf":
        return {
            "compression": "lzf",
            "shuffle": case.get("shuffle", False),
        }
    if kind == "zstd":
        if hdf5plugin is None:
            raise RuntimeError("hdf5plugin is not installed")
        return hdf5plugin.Zstd(clevel=case["level"])
    if kind == "bitshuffle_lz4":
        if hdf5plugin is None:
            raise RuntimeError("hdf5plugin is not installed")
        return hdf5plugin.Bitshuffle(cname="lz4")
    if kind == "blosc2_lz4":
        if hdf5plugin is None:
            raise RuntimeError("hdf5plugin is not installed")
        return hdf5plugin.Blosc2(cname="lz4", clevel=case["level"], filters=hdf5plugin.Blosc2.BITSHUFFLE)
    if kind == "blosc2_zstd":
        if hdf5plugin is None:
            raise RuntimeError("hdf5plugin is not installed")
        return hdf5plugin.Blosc2(cname="zstd", clevel=case["level"], filters=hdf5plugin.Blosc2.BITSHUFFLE)
    if kind == "bitshuffle_zstd":
        if hdf5plugin is None:
            raise RuntimeError("hdf5plugin is not installed")
        return hdf5plugin.Bitshuffle(cname="zstd", clevel=case["level"])
    raise ValueError(f"Unknown case kind: {kind}")


def iter_blocks(shape: tuple[int, int, int], chunk: tuple[int, int, int]):
    npol, nrow, nchan = shape
    cp, cr, cc = chunk
    if cp < npol:
        raise ValueError(f"Chunk polar axis {cp} must cover full polar axis {npol}")
    for row_start in range(0, nrow, cr):
        row_stop = min(row_start + cr, nrow)
        for chan_start in range(0, nchan, cc):
            chan_stop = min(chan_start + cc, nchan)
            yield row_start, row_stop, chan_start, chan_stop


def choose_baseline_chunk(shape: tuple[int, int, int]) -> tuple[int, int, int]:
    npol, nrow, nchan = shape
    return (
        npol,
        min(nrow, 64),
        min(nchan, 1024),
    )


def benchmark_case(source_ds, out_path: Path, case: dict, shape: tuple[int, int, int], raw_bytes: int, baseline_bytes: int) -> dict:
    chunk = tuple(case["chunk"])
    kwargs = resolve_case_kwargs(case)
    write_start = time.perf_counter()
    with h5py.File(out_path, "w") as f:
        ds = f.create_dataset("DATA", shape=shape, dtype=source_ds.dtype, chunks=chunk, **kwargs)
        for row_start, row_stop, chan_start, chan_stop in iter_blocks(shape, chunk):
            ds[:, row_start:row_stop, chan_start:chan_stop] = source_ds[:, row_start:row_stop, chan_start:chan_stop]
    write_seconds = time.perf_counter() - write_start

    file_size = out_path.stat().st_size

    read_start = time.perf_counter()
    with h5py.File(out_path, "r") as f:
        ds = f["DATA"]
        checksum = float(ds[:, :: max(1, shape[1] // 8), :: max(1, shape[2] // 8)].sum(dtype="float64"))
    read_seconds = time.perf_counter() - read_start

    return {
        "name": case["name"],
        "kind": case["kind"],
        "level": case.get("level"),
        "chunk": chunk,
        "write_seconds": write_seconds,
        "read_seconds": read_seconds,
        "file_size_bytes": file_size,
        "raw_data_bytes": raw_bytes,
        "baseline_uncompressed_hdf5_bytes": baseline_bytes,
        "compression_ratio_vs_raw": raw_bytes / file_size if file_size else None,
        "compression_ratio_vs_uncompressed_hdf5": baseline_bytes / file_size if file_size else None,
        "checksum": checksum,
    }


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    selected_cases = select_cases(args.cases)

    with h5py.File(args.source_hdf5, "r") as f:
        source_ds = f["1"]["DATA"]
        shape = list(source_ds.shape)
        if args.max_rows is not None:
            shape[1] = min(shape[1], args.max_rows)
        if args.max_chans is not None:
            shape[2] = min(shape[2], args.max_chans)
        shape = tuple(shape)
        raw_bytes = int(shape[0] * shape[1] * shape[2] * source_ds.dtype.itemsize)

        baseline_path = outdir / "baseline_uncompressed.hdf5"
        baseline_chunk = choose_baseline_chunk(shape)
        write_start = time.perf_counter()
        with h5py.File(baseline_path, "w") as bf:
            ds = bf.create_dataset("DATA", shape=shape, dtype=source_ds.dtype, chunks=baseline_chunk)
            for row_start, row_stop, chan_start, chan_stop in iter_blocks(shape, baseline_chunk):
                ds[:, row_start:row_stop, chan_start:chan_stop] = source_ds[:, row_start:row_stop, chan_start:chan_stop]
        baseline_write_seconds = time.perf_counter() - write_start
        baseline_bytes = baseline_path.stat().st_size

        results = []
        for case in selected_cases:
            out_path = outdir / f"{case['name']}.hdf5"
            result = benchmark_case(source_ds, out_path, case, shape, raw_bytes, baseline_bytes)
            results.append(result)
            print(json.dumps(result, ensure_ascii=False))

    summary = {
        "source_hdf5": args.source_hdf5,
        "shape_used": shape,
        "selected_case_names": [case["name"] for case in selected_cases],
        "raw_data_bytes": raw_bytes,
        "baseline_uncompressed_hdf5_bytes": baseline_bytes,
        "baseline_uncompressed_write_seconds": baseline_write_seconds,
        "results": results,
    }
    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
