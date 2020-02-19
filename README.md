[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3674643.svg)](http://dx.doi.org/10.5281/zenodo.3674643)

# ACC: Auto-Correlogram Calculation in seismology

Extracting P-wave reflections between the free surface and the lithospheric discontinuities to image the subsurface structures.

## Requirements

* Python 3 
* python packages including: 'click', 'commentjson', 'geographiclib', 'matplotlib>=2', 'numpy', 'obspy>=1.0.3', 'pandas', 'setuptools', 'shapely', 'scipy>=0.19.0', 'tqdm'

## Installation

Here I offer two conventional ways to install the package. The first is downloading the code via git clone command.

```
>>> git clone https://github.com/weijias-opensource/acc.git
```

and enter the main directory of the package where the `setup.py` file is, then execute

```
>>> python setup.py install
```
. The second is just simply to execute the command of 

```
>>> pip install seis-acc
```

I strongly suggest you install [Anaconda3](https://docs.anaconda.com/anaconda/install/) first, since 


>Anaconda Distribution is a free, easy-to-install package manager, environment manager, and Python distribution with a collection of 1,500+ open source packages with free community support. Anaconda is platform-agnostic, so you can use it whether you are on Windows, macOS, or Linux.


This allow you to install the acc package using the second way above easily.


## Tutorials

Please go the the example directory and run 

```
>>> sh run.sh
``` 

for a simple example of the Warramungga array data. 

More information can be found at https://acc.readthedocs.io/en/latest/index.html.


## Deployment

The package could be run on all operating systems, including Windows, Mac and Linux. But the package is well-tested on Ubuntu Linux (19.10) with Anaconda3 at now.

## Authors

* **Weijia Sun**

If you have any suggestions to help improve the package, please let me know and I will try to carry them out as soon.

## Contributors

* **B. L. N. Kennett**
* **Huaiyu Yuan**

## Acknowledgments

The author learned to write a flexible, practical and modern software for friendly usage from other packages. More, a small part of code in this package is also reproduced from other projects.

* [rf](https://github.com/trichter/rf)
* [seispy](https://github.com/xumi1993/seispy)
* [obspy](https://github.com/obspy/obspy)

