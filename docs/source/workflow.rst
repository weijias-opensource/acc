Work flow
=========

There are five steps to run the acc package including

0. config the parameter file

1. import raw data of teleseismic events

2. calculate auto-correlograms

3. migrate the P-wave reflectivities to depth domain

4. plot results

0. config parameter
-------------------
Before you start to config parameter, it is suggested to know the principles of the approach. 
Please refer to the following papers.

::

    1. Weijia Sun and B. L. N. Kennett, 2016, Receiver structure from teleseisms: Autocorrelation and cross correlation, Geophys Res Lett, 43, 6234-6242.
    2. Weijia Sun and B. L. N. Kennett, 2017, Mid-lithosphere discontinuities beneath the western and central North China Craton, Geophys Res Lett, 44, 1302-1310. 

Then please read through the configure file carefully. 
Make sure each parameter has been well set. 
For example, the raw data should be in 
**SAC format with the headers of event and station information**.
Then you can set `data` to the path where the data is saved. Wildcard is supported.


1. other steps
--------------

If the parameters are given properly, 
the other four steps can be simply executed via the following commands.


    >>> acc importdata -p conf.json
    >>> acc calevent -p conf.json
    >>> acc migration -p conf.json
    >>> acc plotprofile -p conf.json


The figures will be saved under path_you_set/figures in both pdf and png format.

