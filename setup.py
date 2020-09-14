
from setuptools import setup
from debarcer.src.version import __version__



# Utility function to read the README file.

with open("README.md") as infile:
    content = infile.read().rstrip()

setup(
	name = "debarcer",
	version = __version__,
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
    entry_points={'console_scripts': ['debarcer = debarcer']},
	install_requires = ["numpy>=1.14.2", "pandas>=0.22.0", "pysam>=0.14.1",
                        "matplotlib>=3.1.2", "mistune>=0.8.4", "networkx>=1.11",
                        "pygal>=2.4.0", "scipy>=1.0.1", "seaborn>=0.9.0",
                        "setuptools>=1.1", "pyyaml>=5.1.1"],
    python_requires=">=3.6",
)