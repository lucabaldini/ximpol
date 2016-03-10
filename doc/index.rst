.. ximpol documentation master file, created by
   sphinx-quickstart on Wed Dec  9 11:52:55 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ximpol: an X-ray polarimetry simulation framework
=================================================

ximpol is a simulation framework specifically developed for X-ray polarimetric
applications, based on the `Python <https://www.python.org/>`_ programming
language and the `SciPy <http://www.scipy.org/>`_ stack.
It is not tied to any specific mission or instrument design and is meant
to produce fast and yet realistic observation-simulations, given as basic
inputs:

* an arbitrary source model including morphological, temporal, spectral and
  polarimetric information;
* the response functions of the detector under study, i.e., the effective area,
  the energy dispersion, the point-spread function and the modulation factor.

The format of the response files is OGIP compliant, and the framework has the
capability of producing output files that can be directly fed into the standard
visualization and analysis tools used by the X-ray community, including
XSPEC---which make it a useful tool not only for simulating physical
systems, but also to develop and test end-to-end analysis chains.

If you are wondering what this is all about, you might want to start off by
taking a look at our :ref:`showcase <showcase>`.


Contents:
---------

.. toctree::
   :maxdepth: 2

   showcase
   quick_start
   installation
   source_models
   response_functions
   code_development
   system_tests
   release_notes
   about



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

