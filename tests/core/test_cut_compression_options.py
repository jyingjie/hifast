import unittest
from unittest.mock import patch

from hifast.utils.h5compression import (
    H5CompressionConfig,
    resolve_h5_chunk,
    resolve_h5_dataset_kwargs,
)


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


if __name__ == "__main__":
    unittest.main()
