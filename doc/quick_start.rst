Quick start
===========

The main purpose of this simulation package is to simulate an observation
of a given source model, based on suitable detector response functions. 

The main Monte Carlo simulation script is `ximpol/bin/xpobssim
<https://github.com/lucabaldini/ximpol/blob/master/ximpol/bin/xpobssim.py>`_
and its signature is

.. code-block:: bash

    ximpol/bin/xpobssim.py 
    usage: xpobssim.py [-h] -o OUTPUT_FILE -c CONFIG_FILE [-r IRF_NAME]
                       [-d DURATION] [-t START_TIME] [-n TIME_STEPS]
                       [-s RANDOM_SEED]

Assuming that the current working directory is the ximpol root folder, the
command

.. code-block:: bash

    ximpol/bin/xpobssim.py -c ximpol/srcmodel/config/stationary_point_pl.py -d 100 -o test.fits

should produce an event (FITS) file with a 100 s simulation of a stationary
source with a power-law spectrum (with an index of 2 and normalization of 10) with energy- and time-independent polarization
degree and angle (correctly folded with all the instrument response functions:
effective area, modulation factor, energy dispersion and point-spread function).
The format definition for the event file is in `ximpol/evt/event.py
<https://github.com/lucabaldini/ximpol/blob/master/ximpol/evt/event.py>`_.

You can take a quick look at the output file by typing

.. code-block:: bash

    ximpol/bin/xpevtview.py test.fits

We are already fully equipped for a basic spectral analysis with XSPEC. The
first step is to bin the event file by running the xpbin tool (which creates a `test.pha` file):

.. code-block:: bash

    ximpol/bin/xpbin.py test.fits

Finally we can feed the binned file (along with the corresponding .arf and .rmf response functions) into XSPEC and recover the input parameters of our source.

.. code-block:: bash

    ximpol/bin/xpxspec.py test.pha

Note that xpspec.py is an example of how to use `pyXspec <https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/index.html>`_, unfortunately not all of the XSPEC capabilties have been implemented in pyXspec (for example how to save a plot) so it is left to the user to decide whether to use XSPEC or pyXSPEC for the spectral analysis.

Below is the output from XSPEC on test.pha:

.. code-block:: bash

    Model powerlaw<1> Source No.: 1   Active/On
    Model Model Component  Parameter  Unit     Value
     par  comp
      1    1   powerlaw   PhoIndex            2.00546      +/-  9.41951E-03  
      2    1   powerlaw   norm                10.0265      +/-  7.12876E-02  


    Test statistic : Chi-Squared =         196.87 using 220 PHA bins.
    Reduced chi-squared =        0.90308 for    218 degrees of freedom 
    Null hypothesis probability =   8.447913e-01


.. image:: figures/xspec_screenshot.png
