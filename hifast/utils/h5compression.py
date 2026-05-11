from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class H5CompressionConfig:
    compression: str | int | None = "none"
    compression_level: int | None = None
    chunk_rows: int | None = None
    chunk_chans: int | None = None


H5_COMPRESSION_HELP = (
    'HDF5 compression method. Choices: `none`, `gzip`, `lzf`, '
    '`blosc2_lz4`, `blosc2_zstd`, `bitshuffle_lz4`, `bitshuffle_zstd`, or a number in range(10) as a '
    'backward-compatible alias for `gzip` level.\n'
    'Recommended combinations:\n'
    '  Compatibility/CARTA: `--h5_compression gzip --h5_compression_level 2 '
    '--h5_chunk_rows 128 --h5_chunk_chans 512`\n'
    '  Higher compression: `--h5_compression gzip --h5_compression_level 7 '
    '--h5_chunk_rows 128 --h5_chunk_chans 512`\n'
    '  Faster plugin mode: `--h5_compression blosc2_lz4 '
    '--h5_chunk_rows 128 --h5_chunk_chans 1024`\n'
    '  Balanced plugin mode: `--h5_compression blosc2_zstd '
    '--h5_compression_level 5 --h5_chunk_rows 128 --h5_chunk_chans 1024`\n'
    'Sample benchmark on recent full-frequency 2-pol outputs:\n'
    '  gzip level 2, 128x512: about 1.50x compression, write ~62-80 s, sampled read ~0.03-0.04 s\n'
    '  gzip level 7, 128x512: about 1.51x compression, write ~94-104 s, sampled read ~0.09 s\n'
    '  blosc2_lz4, 128x1024: about 1.47x compression, write ~9 s, sampled read ~0.05 s\n'
    '  blosc2_zstd level 5, 128x1024: about 1.50x compression, write ~15 s, sampled read ~0.04 s\n'
    '  Read timings above come from the benchmark script''s sampled read path rather than full-file sequential reads,\n'
    '  so use them only as a rough comparison between methods.\n'
    'Note: plugin compression needs `hdf5plugin`; when opening plugin-compressed files '
    'in other HDF5 programs, you may need to set `HDF5_PLUGIN_PATH` to the plugin directory. '
    'Prebuilt plugin binaries can be taken from the `hdf5plugin` package: '
    'https://pypi.org/project/hdf5plugin/'
)

H5_COMPRESSION_LEVEL_HELP = (
    'compression level for `gzip`, `blosc2_lz4`, `blosc2_zstd`, or `bitshuffle_zstd`; '
    'default `2` for gzip and `5` for plugin zstd/blosc2 modes'
)

H5_CHUNK_ROWS_HELP = (
    'chunk size along the row axis of output datasets; method-dependent default if omitted'
)

H5_CHUNK_CHANS_HELP = (
    'chunk size along the channel axis of output datasets; method-dependent default if omitted'
)


def load_hdf5plugin():
    try:
        import hdf5plugin
    except ImportError as exc:
        raise ImportError("plugin compression needs `hdf5plugin` to be installed") from exc
    return hdf5plugin


def add_h5_compression_arguments(parser):
    parser.add_argument('--h5_compression', default='none', env_var='HIFAST_H5_COMPRESSION', help=H5_COMPRESSION_HELP)
    parser.add_argument('--h5_compression_level', type=int, env_var='HIFAST_H5_COMPRESSION_LEVEL', help=H5_COMPRESSION_LEVEL_HELP)
    parser.add_argument('--h5_chunk_rows', type=int, env_var='HIFAST_H5_CHUNK_ROWS', help=H5_CHUNK_ROWS_HELP)
    parser.add_argument('--h5_chunk_chans', type=int, env_var='HIFAST_H5_CHUNK_CHANS', help=H5_CHUNK_CHANS_HELP)


def normalize_h5_compression_args(args):
    try:
        legacy_level = int(args.h5_compression)
    except ValueError:
        legacy_level = None
    if legacy_level is not None:
        if not 0 <= legacy_level <= 9:
            raise ValueError('legacy numeric --h5_compression should be in range(10)')
        if args.h5_compression_level is not None:
            raise ValueError('do not set both numeric --h5_compression and --h5_compression_level')
        args.h5_compression = 'gzip'
        args.h5_compression_level = legacy_level
    if args.h5_chunk_rows is not None and args.h5_chunk_rows <= 0:
        raise ValueError('--h5_chunk_rows should be a positive integer')
    if args.h5_chunk_chans is not None and args.h5_chunk_chans <= 0:
        raise ValueError('--h5_chunk_chans should be a positive integer')
    return H5CompressionConfig(
        compression=args.h5_compression,
        compression_level=args.h5_compression_level,
        chunk_rows=args.h5_chunk_rows,
        chunk_chans=args.h5_chunk_chans,
    )


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
    dims = tuple(int(v) for v in shape)
    if len(dims) == 0:
        return None
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
    if len(dims) == 1:
        return (min(dims[0], chunk_chans),)
    if len(dims) == 2:
        return (min(dims[0], chunk_rows), min(dims[1], chunk_chans))
    prefix = dims[:-2]
    return prefix + (min(dims[-2], chunk_rows), min(dims[-1], chunk_chans))


def resolve_h5_dataset_kwargs(shape, config: H5CompressionConfig):
    method = normalize_h5_compression_method(config.compression)
    level = config.compression if isinstance(config.compression, int) else config.compression_level
    chunk = resolve_h5_chunk(
        shape,
        method,
        chunk_rows=config.chunk_rows,
        chunk_chans=config.chunk_chans,
    )

    kwargs = {}
    if chunk is not None:
        kwargs["chunks"] = chunk
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
