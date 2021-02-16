# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys

from gaga2 import __version__ as tool_version

here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    description = long_description = description.read()

    name="gaga2"
    version = tool_version #[line.strip().split(" ")[-1] for line in open("gaga2/__init__.py") if line.startswith("__version__")][0]

    if sys.version_info.major != 3:
        raise EnvironmentError("""{toolname} is a python module that requires python3, and is not compatible with python2.""".format(toolname=name))

    setup(
        name=name,
        version=version,
        description=description,
        long_description=long_description,
        url="https://github.com/cschu/gaga2",
        author="Christian Schudoma",
        author_email="christian.schudoma@embl.de",
        license="MIT",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Topic :: Scientific Engineering :: Bio/Informatics",
            "License :: OSI Approved :: MIT License",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python :: 3.7"
        ],
        zip_safe=False,
        keywords="16S amplicon sequencing processing DADA2 pipeline",
        packages=find_packages(exclude=["test"]),
        #scripts=["util/gff_indexer.py"],
        install_requires=[
        ],
        entry_points={
            "console_scripts": [
                "check_readsets=gaga2.check_readsets:main",
                "gather_fastq_files=gaga2.gather_fastq_files:main",
                "trim_params=gaga2.trim_params:main"
            ],
        },
        scripts=[
            "nextflow/gaga2.nf",
            "R_scripts/dada2_analysis.R",
            "R_scripts/dada2_preprocess.R",
        ],
        package_data={
            "gaga2.etc": ["*.yml"]
        },
        include_package_data=True,
        data_files=[],
    )
