import unittest
import tempfile
from unittest.mock import patch

import numpy as np
import h5py

from hifast.utils.h5compression import (
    H5CompressionConfig,
    add_h5_compression_arguments,
    normalize_h5_compression_args,
    resolve_h5_chunk,
    resolve_h5_dataset_kwargs,
)
from hifast.utils.io import save_specs_hdf5
from hifast.utils.io import BaseIO
from hifast.utils.io import ArgumentParser


class _FakePlugin:
    @staticmethod
    def Bitshuffle(**kwargs):
        return {
            "fake_filter": "bitshuffle",
            "plugin_kwargs": kwargs,
        }


class TestCutCompressionOptions(unittest.TestCase):
    def test_resolve_h5_chunk_defaults(self):
        shape = (2, 2048, 65536)

        self.assertEqual(resolve_h5_chunk(shape, "gzip"), (2, 128, 512))
        self.assertEqual(resolve_h5_chunk(shape, "lzf"), (2, 64, 1024))
        self.assertEqual(resolve_h5_chunk(shape, "blosc2_lz4"), (2, 128, 1024))
        self.assertEqual(resolve_h5_chunk(shape, "blosc2_zstd"), (2, 128, 1024))
        self.assertEqual(resolve_h5_chunk(shape, "bitshuffle_lz4"), (2, 128, 1024))
        self.assertEqual(resolve_h5_chunk(shape, "bitshuffle_zstd"), (2, 128, 1024))

    def test_resolve_h5_chunk_clips_to_dataset_shape(self):
        shape = (2, 100, 400)

        self.assertEqual(resolve_h5_chunk(shape, "gzip"), (2, 100, 400))
        self.assertEqual(resolve_h5_chunk(shape, "bitshuffle_lz4"), (2, 100, 400))

    def test_resolve_h5_chunk_for_2d_and_1d(self):
        self.assertEqual(resolve_h5_chunk((2048, 65536), "gzip"), (128, 512))
        self.assertEqual(resolve_h5_chunk((65536,), "gzip"), (512,))

    def test_resolve_gzip_kwargs_defaults(self):
        kwargs = resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="gzip"))

        self.assertEqual(kwargs["chunks"], (2, 128, 512))
        self.assertEqual(kwargs["compression"], "gzip")
        self.assertEqual(kwargs["compression_opts"], 2)
        self.assertTrue(kwargs["shuffle"])

    def test_resolve_legacy_numeric_gzip_level(self):
        kwargs = resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression=3))

        self.assertEqual(kwargs["chunks"], (2, 128, 512))
        self.assertEqual(kwargs["compression"], "gzip")
        self.assertEqual(kwargs["compression_opts"], 3)
        self.assertTrue(kwargs["shuffle"])

    def test_resolve_bitshuffle_lz4_kwargs(self):
        with patch("hifast.utils.h5compression.load_hdf5plugin", return_value=_FakePlugin):
            kwargs = resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="bitshuffle_lz4"))

        self.assertEqual(kwargs["chunks"], (2, 128, 1024))
        self.assertEqual(kwargs["fake_filter"], "bitshuffle")
        self.assertEqual(kwargs["plugin_kwargs"], {"cname": "lz4"})

    def test_resolve_bitshuffle_zstd_kwargs(self):
        with patch("hifast.utils.h5compression.load_hdf5plugin", return_value=_FakePlugin):
            kwargs = resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="bitshuffle_zstd"))

        self.assertEqual(kwargs["chunks"], (2, 128, 1024))
        self.assertEqual(kwargs["fake_filter"], "bitshuffle")
        self.assertEqual(kwargs["plugin_kwargs"], {"cname": "zstd", "clevel": 5})

    def test_resolve_blosc2_lz4_kwargs(self):
        class _FakeBlosc2Plugin(_FakePlugin):
            class Blosc2:
                BITSHUFFLE = 2

                def __new__(cls, **kwargs):
                    return {
                        "fake_filter": "blosc2",
                        "plugin_kwargs": kwargs,
                    }

        with patch("hifast.utils.h5compression.load_hdf5plugin", return_value=_FakeBlosc2Plugin):
            kwargs = resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="blosc2_lz4"))

        self.assertEqual(kwargs["chunks"], (2, 128, 1024))
        self.assertEqual(kwargs["fake_filter"], "blosc2")
        self.assertEqual(kwargs["plugin_kwargs"], {"cname": "lz4", "clevel": 5, "filters": 2})

    def test_resolve_blosc2_zstd_kwargs(self):
        class _FakeBlosc2Plugin(_FakePlugin):
            class Blosc2:
                BITSHUFFLE = 2

                def __new__(cls, **kwargs):
                    return {
                        "fake_filter": "blosc2",
                        "plugin_kwargs": kwargs,
                    }

        with patch("hifast.utils.h5compression.load_hdf5plugin", return_value=_FakeBlosc2Plugin):
            kwargs = resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="blosc2_zstd"))

        self.assertEqual(kwargs["chunks"], (2, 128, 1024))
        self.assertEqual(kwargs["fake_filter"], "blosc2")
        self.assertEqual(kwargs["plugin_kwargs"], {"cname": "zstd", "clevel": 5, "filters": 2})

    def test_invalid_none_level(self):
        with self.assertRaisesRegex(ValueError, "`none` does not use h5_compression_level"):
            resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="none", compression_level=1))

    def test_invalid_lzf_level(self):
        with self.assertRaisesRegex(ValueError, "`lzf` does not use h5_compression_level"):
            resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="lzf", compression_level=1))

    def test_invalid_bitshuffle_lz4_level(self):
        with patch("hifast.utils.h5compression.load_hdf5plugin", return_value=_FakePlugin):
            with self.assertRaisesRegex(ValueError, "`bitshuffle_lz4` does not use h5_compression_level"):
                resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="bitshuffle_lz4", compression_level=1))

    def test_invalid_gzip_level(self):
        with self.assertRaisesRegex(ValueError, "gzip compression level should be in range\\(10\\)"):
            resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="gzip", compression_level=10))

    def test_invalid_blosc2_lz4_level(self):
        class _FakeBlosc2Plugin(_FakePlugin):
            class Blosc2:
                BITSHUFFLE = 2

                def __new__(cls, **kwargs):
                    return kwargs

        with patch("hifast.utils.h5compression.load_hdf5plugin", return_value=_FakeBlosc2Plugin):
            with self.assertRaisesRegex(ValueError, "blosc2_lz4 compression level should be in \\[0, 9\\]"):
                resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="blosc2_lz4", compression_level=10))

    def test_invalid_blosc2_zstd_level(self):
        class _FakeBlosc2Plugin(_FakePlugin):
            class Blosc2:
                BITSHUFFLE = 2

                def __new__(cls, **kwargs):
                    return kwargs

        with patch("hifast.utils.h5compression.load_hdf5plugin", return_value=_FakeBlosc2Plugin):
            with self.assertRaisesRegex(ValueError, "blosc2_zstd compression level should be in \\[0, 9\\]"):
                resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="blosc2_zstd", compression_level=10))

    def test_invalid_bitshuffle_zstd_level(self):
        with patch("hifast.utils.h5compression.load_hdf5plugin", return_value=_FakePlugin):
            with self.assertRaisesRegex(ValueError, "bitshuffle_zstd compression level should be in \\[1, 22\\]"):
                resolve_h5_dataset_kwargs((2, 2048, 65536), H5CompressionConfig(compression="bitshuffle_zstd", compression_level=0))

    def test_normalize_h5_compression_args_legacy_numeric(self):
        class _Args:
            h5_compression = "3"
            h5_compression_level = None
            h5_chunk_rows = 64
            h5_chunk_chans = 256

        config = normalize_h5_compression_args(_Args)
        self.assertEqual(config, H5CompressionConfig(compression="gzip", compression_level=3, chunk_rows=64, chunk_chans=256))

    def test_add_h5_compression_arguments_sets_env_vars(self):
        parser = ArgumentParser(allow_abbrev=False)
        add_h5_compression_arguments(parser)
        actions = {action.dest: action for action in parser._actions}
        self.assertEqual(actions["h5_compression"].env_var, "HIFAST_H5_COMPRESSION")
        self.assertEqual(actions["h5_compression_level"].env_var, "HIFAST_H5_COMPRESSION_LEVEL")
        self.assertEqual(actions["h5_chunk_rows"].env_var, "HIFAST_H5_CHUNK_ROWS")
        self.assertEqual(actions["h5_chunk_chans"].env_var, "HIFAST_H5_CHUNK_CHANS")

    def test_save_specs_hdf5_applies_gzip_compression(self):
        dict_in = {
            "Header": {"test": "yes"},
            "Ta": np.zeros((2, 32, 64), dtype=np.float32),
            "freq": np.linspace(1000, 1001, 64),
            "mjd": np.linspace(60000, 60001, 32),
        }
        config = H5CompressionConfig(compression="gzip", compression_level=2)
        with tempfile.NamedTemporaryFile(suffix=".hdf5") as tmp:
            save_specs_hdf5(tmp.name, dict_in, h5_compression_config=config)
            with h5py.File(tmp.name, "r") as f:
                self.assertEqual(f["S"]["Ta"].compression, "gzip")
                self.assertEqual(f["S"]["Ta"].compression_opts, 2)
                self.assertEqual(f["S"]["Ta"].chunks, (2, 32, 64))
                self.assertEqual(f["S"]["freq"].compression, "gzip")
                self.assertEqual(f["S"]["mjd"].compression, "gzip")

    def test_baseio_builds_default_h5_compression_config(self):
        class _Args:
            pass

        io = BaseIO.__new__(BaseIO)
        io.args = _Args()
        self.assertEqual(io._get_h5_compression_config(), H5CompressionConfig())

    def test_baseio_builds_h5_compression_config_from_args(self):
        class _Args:
            h5_compression = "blosc2_zstd"
            h5_compression_level = 5
            h5_chunk_rows = 128
            h5_chunk_chans = 1024

        io = BaseIO.__new__(BaseIO)
        io.args = _Args()
        self.assertEqual(
            io._get_h5_compression_config(),
            H5CompressionConfig(
                compression="blosc2_zstd",
                compression_level=5,
                chunk_rows=128,
                chunk_chans=1024,
            ),
        )


if __name__ == "__main__":
    unittest.main()
