# Copyright 2013-2018 Tom Eulenfeld, MIT license
"""
rf Documentation
================

rf is a Python framework for receiver function analysis.
Read and write support of necessary metadata is provided for
SAC, SeismicHandler and HDF5 waveform files.
Data is handled by the ``RFStream`` class which inherits a lot of useful
methods from its ObsPy ancestor ``Stream``,
but also has some unique methods necessary for receiver function calculation.


Method
------

The receiver function method is a popular technique to investigate crustal and

| *Left*: In a two-layer-model part of the incoming P-wave is converted to a
    S-wave at the layer boundary. Major multiples are Pppp, Ppps and Ppss.
| *Right*: Synthetic receiver function of Q component in a two-layer-model.

Installation
------------

Dependencies of rf are

    * ObsPy_ and some of its dependencies,
    * cartopy, geographiclib, shapely,
    * toeplitz_ (time domain deconvolution), tqdm,
    * obspyh5_ for hdf5 file support (optional).

After the installation of Obspy rf can be installed with ::

    pip install rf


Using the Python module
-----------------------

The main functionality is provided by the class `.RFStream`
which is derived from ObsPy's `~obspy.core.stream.Stream` class.


=================  =========  ======
stats              SH/Q       SAC
=================  =========  ======
station_latitude   COMMENT    stla
station_longitude  COMMENT    stlo
station_elevation  COMMENT    stel
event_latitude     LAT        evla
event_longitude    LON        evlo
event_depth        DEPTH      evdp
event_magnitude    MAGNITUDE  mag
event_time         ORIGIN     o
onset              P-ONSET    a
type               COMMENT    kuser0
phase              COMMENT    kuser1
moveout            COMMENT    kuser2
distance           DISTANCE   gcarc
back_azimuth       AZIMUTH    baz
inclination        INCI       user0
slowness           SLOWNESS   user1
pp_latitude        COMMENT    user2
pp_longitude       COMMENT    user3
pp_depth           COMMENT    user4
box_pos            COMMENT    user5
box_length         COMMENT    user6
=================  =========  ======

.. note::
    Q-file header COMMENT is used for storing some information, because
    the Q format has a shortage of predefined headers.


Miscellaneous
-------------

.. _this: http://www.diss.fu-berlin.de/diss/servlets/MCRFileNodeServlet/FUDISS_derivate_000000014929/dissertation_richter.pdf
.. _ObsPy: http://www.obspy.org/
.. _pip: http://www.pip-installer.org/
.. _obspyh5: https://github.com/trichter/obspyh5/
.. _toeplitz: https://github.com/trichter/toeplitz/
.. _GitHub: https://github.com/trichter/rf/
.. |buildstatus| image:: https://api.travis-ci.org/trichter/rf.png?
    branch=master
   :target: https://travis-ci.org/trichter/rf
"""

__version__ = '0.1.1-dev'

# from rf.profile import get_profile_boxes
# from rf.rfstream import read_rf, RFStream, rfstats
# from rf.util import iter_event_data, IterMultipleComponents

# if 'dev' not in __version__:  # get image for correct version from travis-ci
    # __doc__ = __doc__.replace('branch=master', 'branch=v%s' % __version__)
