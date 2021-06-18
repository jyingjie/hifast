import setuptools
import glob
import versioneer
# with open("README.md", "r") as fh:
#     long_description = fh.read()

sofa_lib = setuptools.Extension("_sofa_c",
                       glob.glob('./src/*.c'),
                       depends=["./src/sofa.h", "./src/sofam.h"],
                       include_dirs=["./src"])

setuptools.setup(
    name="hifast", # Replace with your own username
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Yingjie Jing etc.",
    author_email="2012jyj@gmail.com",
    description="A Python module containing functions to process fast data",
#     long_description=long_description,
#     long_description_content_type="text/markdown",
#     url="",
    packages=['hifast','hifast/sofa'],
    package_data={
        "hifast": ["data/*.txt", "data/*.json"],
    },
    ext_package='hifast/sofa',
    ext_modules = [sofa_lib],
    scripts = glob.glob('scripts/*.sh'),
    install_requires=['numpy>=1.12','matplotlib','scipy','h5py','pandas>=1.0','xlrd','astropy','PyAstronomy',
                     'tqdm'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    python_requires='>=3.6',
    zip_safe=False,
)
