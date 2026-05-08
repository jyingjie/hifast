from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class H5CompressionConfig:
    compression: str | int | None = "none"
    compression_level: int | None = None
    chunk_rows: int | None = None
    chunk_chans: int | None = None


def load_hdf5plugin():
    try:
        import hdf5plugin
    except ImportError as exc:
        raise ImportError("`bitshuffle_*` compression needs `hdf5plugin` to be installed") from exc
    return hdf5plugin


def normalize_h5_compression_method(h5_compression):
    if h5_compression is None:
        return "none"
    if isinstance(h5_compression, int):
        if not 0 <= h5_compression <= 9:
            raise ValueError("numeric h5_compression should be in range(10)")
        return "gzip"
    method = str(h5_compression).strip().lower()
    if method == "":
        method = "none"
    valid = {
        "none",
        "gzip",
        "lzf",
        "bitshuffle_lz4",
        "bitshuffle_zstd",
        "blosc2_lz4",
        "blosc2_zstd",
    }
    if method not in valid:
        raise ValueError(f"Unsupported h5 compression method: {h5_compression}")
    return method


def resolve_h5_chunk(shape, method, chunk_rows=None, chunk_chans=None):
    npolar, n_rows, n_chans = [int(v) for v in shape]
    if method == "gzip":
        default_rows, default_chans = 128, 512
    elif method in {"bitshuffle_lz4", "bitshuffle_zstd", "blosc2_lz4", "blosc2_zstd"}:
        default_rows, default_chans = 128, 1024
    else:
        default_rows, default_chans = 64, 1024
    chunk_rows = default_rows if chunk_rows is None else int(chunk_rows)
    chunk_chans = default_chans if chunk_chans is None else int(chunk_chans)
    if chunk_rows <= 0 or chunk_chans <= 0:
        raise ValueError("h5 chunk sizes should be positive integers")
    return (npolar, min(n_rows, chunk_rows), min(n_chans, chunk_chans))


def resolve_h5_dataset_kwargs(shape, config: H5CompressionConfig):
    method = normalize_h5_compression_method(config.compression)
    level = config.compression if isinstance(config.compression, int) else config.compression_level
    chunk = resolve_h5_chunk(
        shape,
        method,
        chunk_rows=config.chunk_rows,
        chunk_chans=config.chunk_chans,
    )

    kwargs = {"chunks": chunk}
    if method == "none":
        if level is not None:
            raise ValueError("`none` does not use h5_compression_level")
        return kwargs
    if method == "gzip":
        if level is None:
            level = 2
        level = int(level)
        if not 0 <= level <= 9:
            raise ValueError("gzip compression level should be in range(10)")
        kwargs.update({
            "compression": "gzip",
            "compression_opts": level,
            "shuffle": True,
        })
        return kwargs
    if method == "lzf":
        if level is not None:
            raise ValueError("`lzf` does not use h5_compression_level")
        kwargs.update({
            "compression": "lzf",
            "shuffle": True,
        })
        return kwargs

    hdf5plugin = load_hdf5plugin()
    if method == "bitshuffle_lz4":
        if level is not None:
            raise ValueError("`bitshuffle_lz4` does not use h5_compression_level")
        kwargs.update(hdf5plugin.Bitshuffle(cname="lz4"))
        return kwargs
    if method == "bitshuffle_zstd":
        if level is None:
            level = 5
        level = int(level)
        if not 1 <= level <= 22:
            raise ValueError("bitshuffle_zstd compression level should be in [1, 22]")
        kwargs.update(hdf5plugin.Bitshuffle(cname="zstd", clevel=level))
        return kwargs
    if method == "blosc2_lz4":
        if level is None:
            level = 5
        level = int(level)
        if not 0 <= level <= 9:
            raise ValueError("blosc2_lz4 compression level should be in [0, 9]")
        kwargs.update(hdf5plugin.Blosc2(cname="lz4", clevel=level, filters=hdf5plugin.Blosc2.BITSHUFFLE))
        return kwargs
    if method == "blosc2_zstd":
        if level is None:
            level = 5
        level = int(level)
        if not 0 <= level <= 9:
            raise ValueError("blosc2_zstd compression level should be in [0, 9]")
        kwargs.update(hdf5plugin.Blosc2(cname="zstd", clevel=level, filters=hdf5plugin.Blosc2.BITSHUFFLE))
        return kwargs
    raise ValueError(f"Unsupported h5 compression method: {method}")
