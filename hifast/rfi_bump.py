"""HiFAST command interface for quasi-periodic narrow bump RFI."""

from __future__ import annotations

__all__ = ["IO", "parser"]

from .utils.io import *

import json
import tempfile
from copy import deepcopy
from pathlib import Path

from .core.rfi_bump import ALGORITHM, apply_masks, detect_masks


def create_parser():
    parser = ArgumentParser(
        prog=(
            f"python -m hifast."
            f"{os.path.basename(sys.argv[0])[:-3]}"
        ),
        formatter_class=formatter_class,
        allow_abbrev=False,
        description=(
            "Experimental pre-Doppler mask for the approximately "
            "1 MHz quasi-periodic narrow bump RFI."
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
            "intermediate-product directory; setting it also "
            "enables keep_products"
        ),
    )
    return parser


parser = create_parser()


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
        output = Path(self.fpath_out)
        return output.with_suffix("").with_name(
            f"{output.stem}_products"
        )

    def _locate_masks(self, workdir):
        start, stop = self._time_range()
        polarizations = (
            ["XX", "YY"]
            if self.args.polar == "both"
            else [self.args.polar]
        )
        return detect_masks(
            fpath=self.args.fpath,
            beam_number=self.nB,
            date=self.date,
            polarizations=polarizations,
            start=start,
            stop=stop,
            workdir=workdir,
            source_frequency_mhz=self.freq,
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
