
from setuptools import setup

# Utility function to read the README file.

with open("README.rst") as infile:
    content = infile.read().rstrip()


version = "2.0.3"

setup(
	name = "debarcer",
	version = version,
	author = "Richard Jovelin",
	author_email = "richard.jovelin@oicr.on.ca",
	description = ("A package for de-barcoding and error correction of sequencing"
											" data containing molecular barcodes."),
	license = "MIT License",
	keywords = "computational genomics",
	url = "https://github.com/oicr-gsi/debarcer",
	packages = ['debarcer', 'tests', 'debarcer.src'],
	long_description = content,
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
	install_requires=['argparse', 'configparser', 'pysam', 'pandas',
                   'matplotlib', 'numpy', 'scipy'],
	python_requires='>=3.6',
)