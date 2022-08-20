import setuptools
import glob
import versioneer

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
    
# with open("README.md", "r") as fh:
#     long_description = fh.read()

setuptools.setup(
    name="hifast", # Replace with your own username
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Yingjie Jing etc.",
    author_email="2012jyj@gmail.com",
    description="https://hifast.readthedocs.io",
#     long_description=long_description,
#     long_description_content_type="text/markdown",
#     url="",
    packages=['hifast', 'hifast.utils', 'hifast.core', 'hifast.ripple', 'hifast.interaction'],
    package_data={
        "hifast.core": ["data/*.txt", "data/*.json"],
    },
    scripts = glob.glob('scripts/*.sh'),
    install_requires=['numpy>=1.12',
                      'matplotlib',
                      'scipy',
                      'h5py',
                      'pandas>=1.0',
                      'openpyxl',
                      'astropy',
                      'configargparse',
                      'tqdm',
                      'threadpoolctl',
                     ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    python_requires='>=3.6',
    zip_safe=False,
)
