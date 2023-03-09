import os
import sys
import time
import subprocess
import argparse
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

parser = argparse.ArgumentParser(prog=f"python -m hifast", formatter_class=argparse.ArgumentDefaultsHelpFormatter, allow_abbrev=False, add_help=True,
                        description='', usage="python -m hifast subcomand fpath1 [fpath2 fpath3 ...] ... [-p P] [Args with '-'] ... [Args with '--'] ...")
parser.add_argument('--version', action='store_true',
                   help='display version')
subparsers = parser.add_subparsers(help='', dest="subcomand", description='')

for subcomand in subcomands_type1:
### time-consuming
#     try:
#         parser_ = import_module(f'.{subcomand}', 'hifast').parser
#         help_s = parser_.format_help().split("positional arguments")
#         epilog = "optional arguments" + help_s[1].split("optional arguments")[1]
#     except Exception as e:
#            print(e)
#
    epilog = f"RUN ``python -m hifast.{subcomand} -h `` to get more info about Args with '-' and '--'"
    usage = f"python -m hifast {subcomand} fpath1 [fpath2 fpath3 ...] ... [-p P] [Args with '-'] ... [Args with '--'] ..."
    _help = f"Type1: RUN ``python -m hifast.{subcomand}`` with all the ``Args with '-'`` and ``Args with '--'`` for file ``fpath1 fpath2 fpath3 ...``"
    parser_a = subparsers.add_parser(subcomand, formatter_class=argparse.RawDescriptionHelpFormatter,
                 help=_help,
                 description=_help,
                 usage=usage,
                 epilog=epilog,
                        )
    parser_a.add_argument('fpath', nargs="+", help='input files path: fpath1 [fpath2 fpath3 ...]')
    parser_a.add_argument('-p', type=int, default=1,
                    help='max processes at a time')
for subcomand in subcomands_type2:
### time-consuming
#     try:
#         parser_ = import_module(f'.{subcomand}', 'hifast').parser
#         help_s = parser_.format_help().split("positional arguments")
#         epilog = "positional arguments" + help_s[1]
#     except:
    epilog = f"RUN ``python -m hifast.{subcomand} -h `` to get more info about Args with '-' and '--'"
    usage = f"python -m hifast {subcomand} fpath1 [fpath2 fpath3 ...] ... [Args with '-'] ... [Args with '--'] ..."
    _help = f"Type2: RUN ``python -m hifast.{subcomand} fpath1 [fpath2 fpath3 ...] ... [Args with '-'] ... [Args with '--'] ..."
    parser_a = subparsers.add_parser(subcomand, formatter_class=argparse.RawDescriptionHelpFormatter,
                 help=_help,
                 description=_help,
                 usage=usage,
                 epilog=epilog,
                                )

## RUN
args, remain = parser.parse_known_args()

if args.version:
    from .__init__ import __version__
    print(__version__)
    sys.exit(0)

# if '--outdir' in remain:
#     ind = remain.index('--outdir') + 1
#     remain[ind] = f"'{remain[ind]}'"

for i in range(len(remain)):
    if '%(' in remain[i] and not remain[i].startswith("'") and not remain[i].startswith("\""):
        remain[i] = f"'{remain[i]}'"

if args.subcomand is None:
    command = f"{sys.executable} -m hifast -h"
elif args.subcomand in subcomands_type1:
    command = f"echo \"{' '.join(args.fpath)}\" | xargs -n 1 -P {args.p} {sys.executable} -m hifast.{args.subcomand} {' '.join(remain)}"
elif args.subcomand in subcomands_type2:
    command = f"{sys.executable} -m hifast.{args.subcomand} {' '.join(remain)}"

p = subprocess.Popen(
        command,
        stderr=sys.stderr, stdout=sys.stdout,
        shell=True
          )
try:
    p.wait()
except KeyboardInterrupt:
    p.kill()
