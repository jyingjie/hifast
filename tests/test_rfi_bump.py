import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from hifast.core.rfi_bump import (
    ALGORITHM,
    CONFIG_PATH,
    apply_masks,
    nearest_frequency_indices,
)
from hifast.rfi_bump import (
    IO,
    parser,
)


def make_input(path):
    flux = np.arange(2 * 5 * 8, dtype=np.float32).reshape(
        2, 5, 8
    )
    with h5py.File(path, "w") as handle:
        header = handle.create_group("Header")
        header.attrs["existing"] = "yes"
        spectra = handle.create_group("S")
        spectra.create_dataset("flux", data=flux)
        spectra.create_dataset(
            "freq",
            data=np.linspace(1420.0, 1425.0, 8),
        )
        spectra.create_dataset(
            "mjd",
            data=np.arange(5, dtype=float),
        )
        spectra.create_dataset(
            "is_rfi",
            data=np.zeros((5, 8), dtype=bool),
        )
    return flux


def result_row(polar):
    return {
        "polar": polar,
        "status": "masked",
        "selected_stage": "primary",
        "reason": "accepted",
        "masked_tooth_count": 3,
        "masked_fraction_target_band": 0.05,
        "masked_sample_count": 2,
        "new_nan_count": 0,
    }


def test_packaged_frozen_configuration_exists():
    assert CONFIG_PATH.is_file()
    assert (
        CONFIG_PATH.parent.parent
        / "run_hierarchical_mask_validation.py"
    ).is_file()


def test_parser_uses_standard_common_and_compression_arguments():
    destinations = {
        action.dest
        for action in parser._actions
    }
    assert {
        "fpath",
        "outdir",
        "force",
        "my_config",
        "h5_compression",
    } <= destinations


def test_output_step_name_uses_expected_separator():
    instance = object.__new__(IO)
    assert instance._get_fpart() == "-rfi_bump"


def test_nearest_frequency_indices_requires_contiguous_match():
    source = np.asarray([1420.0, 1420.1, 1420.2, 1420.3])
    target = np.asarray([1420.1, 1420.2])
    assert np.array_equal(
        nearest_frequency_indices(source, target),
        np.asarray([1, 2]),
    )
    with pytest.raises(ValueError, match="more than 1 Hz"):
        nearest_frequency_indices(
            source,
            np.asarray([1420.15]),
        )


def test_apply_masks_copies_input_and_sets_only_selected_values():
    spectra = np.arange(5 * 8 * 2, dtype=float).reshape(
        5, 8, 2
    )
    original = spectra.copy()
    mask = np.asarray([[True, False], [False, True]])
    results = [result_row("XX")]
    output = apply_masks(
        spectra,
        {
            0: (
                np.asarray([1, 2]),
                np.asarray([3, 4]),
                mask,
            )
        },
        results,
    )
    assert np.array_equal(spectra, original)
    expected_nan = np.zeros_like(output, dtype=bool)
    expected_nan[1, 3, 0] = True
    expected_nan[2, 4, 0] = True
    assert np.array_equal(np.isnan(output), expected_nan)
    assert results[0]["new_nan_count"] == 2


def test_standard_io_writes_new_file_and_preserves_input(
    tmp_path,
    monkeypatch,
):
    source = (
        tmp_path
        / "project-M05_W-20231209-specs_T-flux.hdf5"
    )
    original = make_input(source)
    mask = np.asarray([[True, False], [False, True]])
    results = [result_row("XX")]

    def fake_locate(self, workdir):
        return (
            {
                0: (
                    np.asarray([1, 2]),
                    np.asarray([3, 4]),
                    mask,
                )
            },
            results,
        )

    monkeypatch.setattr(IO, "_locate_masks", fake_locate)
    arguments = parser.parse_args(
        [
            str(source),
            "--outdir",
            str(tmp_path),
        ]
    )
    io = IO(arguments)
    io()

    output = (
        tmp_path
        / "project-M05_W-20231209-specs_T-flux-rfi_bump.hdf5"
    )
    assert output.is_file()
    with h5py.File(source, "r") as handle:
        assert np.array_equal(handle["S/flux"][()], original)
        assert "rfi_bump_algorithm" not in handle["Header"].attrs
    with h5py.File(output, "r") as handle:
        actual = handle["S/flux"][()]
        assert np.isnan(actual[0, 1, 3])
        assert np.isnan(actual[0, 2, 4])
        assert not np.any(handle["S/is_rfi"][()])
        assert handle["Header"].attrs["existing"] == "yes"
        assert (
            handle["Header"].attrs["rfi_bump_algorithm"]
            == ALGORITHM
        )
        recorded = json.loads(
            handle["Header"].attrs["rfi_bump_results"]
        )
        assert recorded[0]["new_nan_count"] == 2
        assert any(
            key.startswith("HISTORY-")
            for key in handle["Header"].attrs
        )
