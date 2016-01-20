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
source with a power-law spectrum with energy- and time-independent polarization
degree and angle (correctly folded with all the instrument response functions:
effective area, modulation factor, energy dispersion and point-spread function).
The format definition for the event file is in `ximpol/evt/event.py
<https://github.com/lucabaldini/ximpol/blob/master/ximpol/evt/event.py>`_.

You can take a quick look at the output file by typing

.. code-block:: bash

    ximpol/bin/xpevtview.py test.fits

We are already fully equipped for a basic spectral analysis with XSPEC. The
first step is to bin the event file

.. code-block:: bash

    ximpol/bin/xpbin.py test.fits

(This will create a corresponding `test.pha` file.) Finally we can feed the
binned file (along with the corresponding .arf and .rmf response functions)
into XSPEC and recover the input parameter of the power law.

.. code-block:: bash

    ximpol/bin/xpxspec.py test.pha

.. image:: figures/xspec_screenshot.png
