import os
import sys
import time
import subprocess
import argparse
import shlex
from glob import glob


# from importlib import import_module

# subcomands_type1 = [
#     'sep',
#     'pos_swi',
#     'radec',
#     'bld',
#     'flux',
#     'sw',
#     'rfi',
#     'multi',
#     'add_radec',
#     'sub_ref',
#     'convert',
#     'downsample',
#      ]

# subcomands_type2 = [
#     
#     'waterfall',
#      ]

subcomands_excl = ['cube', 'cube2', 'waterfall', 'funcs', 'find']
subcomands_type1 = []
for fpath in glob(os.path.dirname(__file__) + '/[a-z,A-Z]*.py'):
    subcomand = os.path.basename(fpath)[:-3]
    if subcomand not in subcomands_excl:
        subcomands_type1 += [subcomand, ]
subcomands_type1.sort()
subcomands_type2 = []



## parser

parser = argparse.ArgumentParser(prog=f"python -m hifast", formatter_class=argparse.ArgumentDefaultsHelpFormatter, allow_abbrev=False, add_help=False,
                        description='', usage="python -m hifast subcomand fpath1 [fpath2 fpath3 ...] ... [-p P] [Args with '-'] ... [Args with '--'] ...")
parser.add_argument('-h', '--help', action='store_true',
                   help='')
parser.add_argument('--version', action='store_true',
                   help='display version')
subparsers = parser.add_subparsers(help='', dest="subcomand", description='')

for subcomand in subcomands_type1:
    parser_a = subparsers.add_parser(subcomand, formatter_class=argparse.RawDescriptionHelpFormatter,
                 add_help=False,
                  )
    parser_a.add_argument('fpath', nargs="*", help='input files path: fpath1 [fpath2 fpath3 ...]')
    parser_a.add_argument('-p', type=int, default=1,
                    help='max processes at a time')
    parser_a.add_argument('-h', '--help', action='store_true',
                    help='')
for subcomand in subcomands_type2:
    parser_a = subparsers.add_parser(subcomand, formatter_class=argparse.RawDescriptionHelpFormatter,
                 add_help=False,
                                )
    parser_a.add_argument('-h', '--help', action='store_true',
                   help='')

## RUN
args, remain = parser.parse_known_args()

if args.version:
    from .__init__ import __version__
    print(__version__)
    sys.exit(0)


if args.subcomand is None:
    print(parser.format_help())
    sys.exit(0)
elif args.help:
    subcomand = args.subcomand
    if subcomand in subcomands_type1:
        _help = f"""
usage: python -m hifast {subcomand} fpath1 [fpath2 fpath3 ...] [-p P] [Arguments of hifast.{subcomand} with '-' or '--' ]

Type1: RUN ``python -m hifast.{subcomand}`` for files in parallel with ``P`` processes.

optional arguments:
-p P        max processes at a time
-h, --help

Below is the help from hifast.{subcomand}:
        """
    elif subcomand in subcomands_type2:
        _help = f"""
Same with hifast.{subcomand}.

Below is the help from hifast.{subcomand}:
        """
    command = f"python -m hifast.{args.subcomand} -h"
    print(_help)
    sys.stdout.flush()
    # help_sub = subprocess.check_output(command, shell=True)
    # help_sub = help_sub.decode('utf-8')
    # print(_help, help_sub)
    # sys.exit(0)
elif args.subcomand in subcomands_type1:
    fpaths = '\n'.join(args.fpath)
    command = f"printf '{fpaths}' | xargs -I {{}} -n 1 -P {args.p} {sys.executable} -m hifast.{args.subcomand} {{}} {shlex.join(remain)}"
elif args.subcomand in subcomands_type2:
    command = f"{sys.executable} -m hifast.{args.subcomand} {shlex.join(remain)}"

p = subprocess.Popen(
        command,
        stderr=sys.stderr, stdout=sys.stdout,
        shell=True
          )
try:
    p.wait()
except KeyboardInterrupt:
    p.kill()
