
import sys
import os
from setuptools import setup

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

version = "2.0.0"

setup(
	name = "debarcer",
	version = version,
	author = "Theodore Bodak",
	author_email = "tabodak@edu.uwaterloo.ca",
	description = ("A package for de-barcoding and error correction of sequencing"
						" data containing molecular barcodes."),
	license = "OSI Approved :: MIT License",
	keywords = "computational genomics",
	url = "https://github.com/oicr-gsi/debarcer",
	packages = ['debarcer', 'tests'],
	long_description = read("README.rst"),
	classifiers = [
	"Development Status :: 3 - Alpha",
	"Intended Audience :: Science/Research",
	"Intended Audience :: Developers",
	"License :: OSI Approved :: MIT License",
	"Programming Language :: Python :: 3",
	"Programming Language :: Python :: 3.6",
	"Topic :: Software Development",
	"Topic :: Scientific/Engineering",
	"Operating System :: POSIX",
	"Operating System :: Unix",
	"Operating System :: MacOS",
	],
	install_requires=['argparse', 'configparser', 'pysam', 'pandas'],
	python_requires='>=3.6',
)
