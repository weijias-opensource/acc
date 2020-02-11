
# ACC: Auto-Correlogram Calculation in seismology

Extracting P-wave reflections between the free surface and the lithospheric discontinuities to image the subsurface structures.

### Installation

Here I provide two conventional ways to install the package. The first is downloading the code via git clone command.

```
>>> git clone https://github.com/weijias-opensource/acc.git
```

and enter the main directroy of the package where the `setup.py` file is, then execute

```
>>> python setup.py install
```
. The second is just simply executing the command of 

```
>>> pip install seis-acc
```


## Tutorials

Please go the the example directory and run 

```
>>> sh run.sh
``` 

for a simple example of the Warramungga array data. More tutorials will be added later.


## Deployment

The package could be running on all operating systems, including Windows, Mac and Linux. But the package is well tested on Ubuntu Linux (19.04) at now.

## Authors

* **Weijia Sun**

If you have any suggestions to help improve the package, please let me know and I will try to implement them as soon.

## Contributors

* **B. L. N. Kennett**

## Acknowledgments

The author learned to write a flexible and practical code for friendly usage from other packages. A small portion of code in this packakge is also reproduced from other projects.

* [rf](https://github.com/trichter/rf)
* [seispy](https://github.com/xumi1993/seispy)
* [obspy](https://github.com/obspy/obspy)

