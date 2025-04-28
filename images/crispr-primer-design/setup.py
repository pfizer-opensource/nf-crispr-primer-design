import re
from setuptools import setup, find_namespace_packages

version = re.search(
    	'^__version__\s*=\s*["\'](.*)["\']',
    	open('src/crispr/primer_design/cli.py').read(),
    	re.M
    	).group(1)

setup(
    name="crispr-primer-design",
    version=version,
    packages=find_namespace_packages(where='src'),
    package_dir={'': 'src'},
    scripts=[],
    entry_points={
        'console_scripts': [
            'CRISPR_primer_design=crispr.primer_design.cli:main',
        ],
    },
    install_requires=[
        'pyyaml>=5.1.2',
        'confuse>=1.7.0,<2',
        'primer3-py>=0.6.1,<1',
    ],
    author="Jason Arroyo",
    author_email="Jason.Arroyo@pfizer.com"
)
