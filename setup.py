import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
DESCRIPTION = "qscgrn"

PKGS = find_packages()
PKG_NAME = "QuantumGRN"
PKG_VERSION = '0.0.1'

MAINTAINER = 'Cristhian Roman'
MAINTAINER_EMAIL = 'cristhianromanvicharra@gmail.com'

PYTHON_REQUIRES = ">=3.7"
URL = "https://github.com/cailab-tamu/QuantumGRN"

LICENSE = "MIT"
CLFS = [
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.7",
]

INSTALL_REQUIRES = [
        "cycler==0.11.0",
        "fonttools==4.33.3",
        "igraph==0.9.11",
        "kiwisolver==1.4.3",
        "matplotlib==3.5.2",
        "numpy==1.23.0",
        "packaging==21.3",
        "pandas==1.4.3",
        "Pillow==9.1.1",
        "pycairo==1.21.0",
        "pyparsing==3.0.9",
        "python-dateutil==2.8.2",
        "pytz==2022.1",
        "six==1.16.0",
        "texttable==1.6.4",
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