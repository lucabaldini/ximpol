Architecture
============

The vast majority of the simulation and data preparation facilities implemented
in ximpol are made available through three main *executables*, as illustred in
the block diagram below:

* ``xpobssim``: given a source model and a set of instrument response
  functions, it produces a photon list correponding to a given observation
  time. The persistent incarnation of a photon list (that we call an
  *event file*) is binary FITS table whose format is detailed here.
* ``xpselect``: it allows to select subsamples of photons in a given
  event file, based on the event energy, direction, time or phase, producing a
  new (smaller) event file.
* ``xpbin``: it allows to bin the data in several different flavors, producing
  counts maps and spectra, light curves, phasograms and *modulation cubes*
  (i.e., histograms of the measured azimuthal distributions in multiple
  energy layers).

.. image:: figures/diagrams/ximpol_diagram.png
   :align: center

All the ximpol tools are fully configurable via command-line and the
relative signatures are details here. In addition to this, ximpol provides
a pipeline facility that allow to combine all the aforementioned
functionalities in a Python script.


Implementation details
======================
