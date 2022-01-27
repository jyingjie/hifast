"""
python -m hifast
"""
import os
from glob import glob
clis = glob(os.path.dirname(__file__)+'/cli*.py')
clis = map(lambda x: '  hifast.' + os.path.basename(x)[:-3], clis)
clis = '\n'.join(clis)

prog = f"""Command-line interface:
{clis}"""
print(prog)
# import argparse
# parser = argparse.ArgumentParser(allow_abbrev=False)
# args = parser.parse_args()
