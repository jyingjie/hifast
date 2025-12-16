import setuptools
import glob
import versioneer
import sys
import os

# Check for loose mode environment variable
is_loose_mode = os.environ.get('HIFAST_LOOSE_DEPENDENCIES') == '1'
if is_loose_mode:
    print("WARNING: HIFAST_LOOSE_DEPENDENCIES is set. Using loose dependency constraints.")

ignore_python_check = os.environ.get('HIFAST_IGNORE_PYTHON_VERSION') == '1'

# Runtime Python version check
# Only enforce strict check if NOT in loose mode and NOT explicitly ignored
if not is_loose_mode and not ignore_python_check:
    if sys.version_info < (3, 9) or sys.version_info >= (3, 10):
        sys.exit("Python >= 3.9 and < 3.10 is required for strict mode. Current version: " + sys.version + 
                 "\nTo use loose dependencies (and bypass this check), set HIFAST_LOOSE_DEPENDENCIES=1" +
                 "\nTo just bypass this check (keep strict deps), set HIFAST_IGNORE_PYTHON_VERSION=1")

import re

# from hifast._version import get_versions
# _re_version = re.compile('^__version__\s*=.*$', re.MULTILINE)
# fname = 'hifast/__init__.py'
# __version__ = get_versions()['version']
# version = f'__version__ = "{__version__}"'

# with open(fname, 'r') as f: code = f.read()
# if _re_version.search(code) is None:
#     code = version + "\n" + code
# else:
#     code = _re_version.sub(version, code)
# with open(fname, 'w') as f:
#     f.write(code)
    

long_description = """
HiFAST is a specialized pipeline developed for the calibration and imaging of neutral atomic hydrogen (HI) data from the Five-hundred-meter Aperture Spherical radio Telescope (FAST).

This package provides a comprehensive toolkit for radio astronomers, facilitating the processing of raw observational data into science-ready data cubes. For in-depth information on the pipeline's algorithms and methodologies, please refer to our [publications](https://hifast.readthedocs.io/en/latest/citations.html).

**Documentation**: https://hifast.readthedocs.io"""

# Define dependency lists
strict_requirements = [
    'numpy~=1.21.0',
    'scipy~=1.7.0',
    'astropy~=4.2.1',

    'matplotlib~=3.4.2',
    'h5py~=3.3.0',
    'openpyxl~=3.0.7',
    'configargparse~=1.5.3',
    'tqdm~=4.61.1',
    'threadpoolctl~=3.1.0',
    'requests~=2.25.1',
    'pyerfa~=2.0.0',
    'Pillow~=8.4.0',
    # Transitive dependencies pinned for stability
    # Removed non-essential deps (certifi, chardet, idna, etc.) to let top-level libs manage them
]

loose_requirements = [
    'numpy>=1.17,<2',
    'scipy',
    'astropy>=4.0',

    'matplotlib',
    'h5py',
    'openpyxl',
    'configargparse',
    'tqdm',
    'threadpoolctl',
    'requests',
    'pyerfa',
    'Pillow',
]

strict_interaction_requirements = [
    'jupyterlab~=3.0.16',
    'ipympl~=0.7.0',
    'ipywidgets~=7.6.3',
    'notebook~=6.4.5',
    'ipykernel~=6.4.1',
    'mpl-interactions~=0.18.1',
]

loose_interaction_requirements = [
    'jupyterlab',
    'ipympl',
    'ipywidgets',
    'notebook',
    'ipykernel',
    'mpl-interactions',
]

# Select requirements
install_requires = loose_requirements if is_loose_mode else strict_requirements
interaction_requires = loose_interaction_requirements if is_loose_mode else strict_interaction_requirements
python_requires = '>=3.6' if is_loose_mode else '>=3.9, <3.10'

setuptools.setup(
    name="hifast", # Replace with your own username
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="HiFAST Developers",
    author_email="2012jyj@gmail.com",
    description="A Python-based pipeline for FAST HI data calibration and imaging",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://hifast.readthedocs.io",
    packages=['hifast',
              'hifast.utils',
              'hifast.core',
              'hifast.ripple',
              'hifast.interaction',
              'hifast.cbr'],
    package_data={
        "hifast.core": ["data/*.txt", "data/*.json"],
        "hifast.cbr": ["data/*.json"],
    },
    scripts = glob.glob('scripts/*.sh'),
    install_requires=install_requires,
    extras_require={
        'interaction': interaction_requires
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
    ],
    python_requires=python_requires,
    # zip_safe=False is required because the code uses __file__ to locate data files (e.g. in hifast/core/flux.py)
    zip_safe=False,
)
