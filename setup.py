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

PYTHON_REQUIRES = ">=3.8"
URL = "https://github.com/cailab-tamu/QuantumGRN"

LICENSE = "MIT"
CLFS = [
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
]

INSTALL_REQUIRES = [
        "certifi~=2022.6.15",
        "cffi~=1.15.1",
        "charset-normalizer~=2.1.0",
        "cryptography~=37.0.2",
        "cycler~=0.11.0",
        "dill~=0.3.5.1",
        "fonttools~=4.33.3",
        "idna~=3.3",
        "igraph~=0.9.11",
        "kiwisolver~=1.4.3",
        "matplotlib~=3.5.2",
        "mpmath~=1.2.1",
        "ntlm-auth~=1.5.0",
        "numpy~=1.23.0",
        "packaging~=21.3",
        "pandas~=1.4.3",
        "pbr~=5.9.0",
        "Pillow~=9.2.0",
        "ply~=3.11",
        "psutil~=5.9.1",
        "pycairo~=1.21.0",
        "pycparser~=2.21",
        "pylatexenc~=2.10",
        "pyparsing~=3.0.9",
        "python-dateutil~=2.8.2",
        "pytz~=2022.1",
        "qiskit~=0.37.0",
        "qiskit-aer~=0.10.4",
        "qiskit-ibmq-provider~=0.19.2",
        "qiskit-terra~=0.21.0",
        "requests~=2.28.1",
        "requests-ntlm~=1.1.0",
        "retworkx~=0.11.0",
        "scipy~=1.8.1",
        "six~=1.16.0",
        "stevedore~=3.5.0",
        "sympy~=1.10.1",
        "texttable~=1.6.4",
        "tweedledum~=1.1.1",
        "urllib3~=1.26.9",
        "websocket-client~=1.3.3",
        "websockets~=10.3"
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
