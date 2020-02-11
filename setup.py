# Copyright 2019 Weijia Sun, MIT license
import os.path
import re

from setuptools import find_packages, setup


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


VERSION = find_version('acc', '__init__.py')
DESCRIPTION = 'Auto-Correlogram Calculation in seismology'
# LONG_DESCRIPTION = (
#     'Please look at the project site for tutorials and information.')
with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()


ENTRY_POINTS = {
    'console_scripts': ['acc=acc.main:run',
                        ]}

REQUIRES = ['cartopy', "click", 'commentjson', 'geographiclib',
            'matplotlib>=2', 'numpy',
            'obspy>=1.0.3', "pandas",
            'setuptools', 'shapely', 'scipy>=0.19.0', 'tqdm'
            ]

# EXTRAS_REQUIRE = {
#     'doc': ['sphinx', 'alabaster'],  # and decorator, obspy
#     'h5': ['obspyh5>=0.3']}

CLASSIFIERS = [
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Topic :: Scientific/Engineering :: Physics'
]


setup(name='seis-acc',
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type="text/markdown",
      url='https://github.com/weijias-opensource/acc',
      author='Weijia SUN',
      author_email='weijia_sun@163.com',
      license='MIT',
      packages=find_packages(),
      package_dir={'acc': 'acc'},
      install_requires=REQUIRES,
      # extras_require=EXTRAS_REQUIRE,
      entry_points=ENTRY_POINTS,
      # please note the entry_points
      # The magic is in the entry_points parameter. Below console_scripts, each line identifies one console script.
      # The first part before the equals sign (=) is the name of the script that should be generated,
      # the second part is the import path followed by a colon (:) with the Click command.
      # entry_points={"console_scripts": ['acc=acc.main:run',],},
      include_package_data=True,
      zip_safe=False,
      classifiers=CLASSIFIERS
      )
