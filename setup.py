import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
DESCRIPTION = "qscgrn"

PKGS = find_packages()
PKG_NAME = "QuantumGRN"
PKG_VERSION = '1.0.1'

MAINTAINER = 'Cristhian Roman'
MAINTAINER_EMAIL = 'cristhianromanvicharra@gmail.com'

PYTHON_REQUIRES = ">=3.9"
URL = "https://github.com/cailab-tamu/QuantumGRN"

LICENSE = "MIT"
CLFS = [
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
]

INSTALL_REQUIRES = [
        "igraph>=0.10",
        "matplotlib>=3.6",
        "numpy>=1.23",
        "pandas>=1.5",
        "pycairo>=1.21",
        "qiskit>=1.0",
        "qiskit-aer>=0.13",
    ]

# This call to setup() does all the work
setup(
    name=PKG_NAME,
    version=PKG_VERSION,
    description=DESCRIPTION,
    long_description=README,
    long_description_content_type="text/markdown",
    url=URL,
    author=MAINTAINER,
    author_email=MAINTAINER_EMAIL,
    license=LICENSE,
    python_requires=PYTHON_REQUIRES,
    classifiers=CLFS,
    packages=PKGS,
    include_package_data=True,
    install_requires=INSTALL_REQUIRES,
    # entry_points={
    #     "console_scripts": [
    #         "scTenifold=scTenifold.__main__:app",
    #     ]
    # },
)
