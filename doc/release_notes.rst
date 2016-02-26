Release notes
=============

* Source model for Cas A added.
* First xpselect implementation (issue #51).


*ximpol (0.28.1) - Thu, 25 Feb 2016 16:48:43 -0800*

* Updated installation instructions (issue #64).
  

*ximpol (0.28.0) - Thu, 25 Feb 2016 15:51:26 -0800*

* Phaseograms implemented in xpbin.py (issue #67).


*ximpol (0.27.0) - Thu, 25 Feb 2016 15:31:53 -0800*

* Work started toward the implementation of periodic sources (issue #43).
* New xEphemeris class in ximpol.srcmodel.roi.py
* New xPeriodicPointSource class in ximpol.srcmodel.roi.py
* Some significant refactoring of the spline and rand classes to allow for
  more flexibility.
* Major change to the source model interface---the energy spectrum and
  polarization degree and angle are now passed to the constructor.
* A whole bunch of obsolete stuff removed from ximpol.srcmodel.spectrum
  (issue #64).
* All configuration files reworked according to the new interfaces.


*ximpol (0.26.0) - Tue, 23 Feb 2016 16:42:27 -0800*

* FILE mode implemented for tbinalg (issue #53).


*ximpol (0.25.0) - Tue, 23 Feb 2016 16:33:27 -0800*

* ebinalgs FILE and EQP implemented (issue #56).


*ximpol (0.24.1) - Tue, 23 Feb 2016 15:55:06 -0800*

* Fixed unit tests.


*ximpol (0.24.0) - Fri, 19 Feb 2016 16:14:36 -0800*

* Vignetting now into the effective area tables (but not used in the
  simulation, yet).


*ximpol (0.23.1) - Thu, 18 Feb 2016 15:03:59 -0800*

* More information added to the IRF primary headers (issue #49).


*ximpol (0.23.0) - Thu, 18 Feb 2016 14:56:15 -0800*

* Major refactoring of ximpol/detector/xipe.py to use the new classes
  (issue #49).
* New optics aeff files provided by Fabio committed (but only the on-axis
  values used for the time being).
* XIPE baseline and goal response functions created (only the effective areas
  differ for the time being).


*ximpol (0.22.4) - Mon, 08 Feb 2016 16:34:11 -0800*

* Fix for issue #59.


*ximpol (0.22.3) - Mon, 08 Feb 2016 16:25:59 -0800*

* Fix for issue #58.


*ximpol (0.22.2) - Mon, 08 Feb 2016 15:51:53 -0800*

* Quick polarization analysis routine in place.
* Bug fix in the new code reading the IRFs.


*ximpol (0.22.1) - Mon, 08 Feb 2016 15:11:38 -0800*

* More refactoring of the binning classes.
* Detector, ROI and IR information propagated from the event to the binned
  files (issue #57).


*ximpol (0.22.0) - Fri, 05 Feb 2016 13:56:10 -0800*

* MCUBE mode implemented in xpbin.py


*ximpol (0.21.2) - Thu, 04 Feb 2016 15:41:41 -0800*

* Source model string formatting improved.
* A few minor changes.


*ximpol (0.21.1) - Thu, 04 Feb 2016 14:28:43 -0800*

* Committed a whole bunch of files left out by mistake.


*ximpol (0.21.0) - Thu, 04 Feb 2016 14:27:20 -0800*

* Major refactoring and revamp of xpevtview.py
* New class for tabulated stationary spectra.
* New configuration file for the SgrB complex.
* Spectral data for the SgrA and SgrB complexes.
* New small utility (xpsrccoords.py) to search for source coordinates.


*ximpol (0.20.0) - Thu, 04 Feb 2016 10:43:26 -0800*

* Gaussian disk spatial template implemented.
* A few srcmodel config files renamed.


*ximpol (0.19.1) - Wed, 03 Feb 2016 16:17:09 -0800*

* Updated documentation.


*ximpol (0.19.0) - Wed, 03 Feb 2016 16:12:42 -0800*

* Uniform disk implemented (issue #54).
* Added command-line option to use the MC Ra/Dec for xpbin.


*ximpol (0.18.0) - Wed, 03 Feb 2016 15:13:52 -0800*

* More work on xpbin.py (closing issues #42 and #52).


*ximpol (0.17.0) - Tue, 02 Feb 2016 15:41:14 -0800*

* Major refactoring of xpbin.py (issue #42).
* Minimum and maximum valid times added to the model components.
* Configuration file for a GRB added.


*ximpol (0.16.1) - Tue, 26 Jan 2016 18:49:19 -0800*

* Minor refactoring of the ximpol.core.fitsio module.
  

*ximpol (0.16.0) - Tue, 26 Jan 2016 18:40:11 -0800*

* Module ximpol.core.fitsio added (issue #49).
* ximpol.evt.event refactored to use the new ximpol.core.fitsio module.
* GTI list in the output event file (issue #24)
* ROI source table in the output event file (issue #45).
* IRF name added in the output event file header (issue #24).
* ROI information added in the output event file header (issue #48).


*ximpol (0.15.2) - Mon, 25 Jan 2016 18:04:33 -0800*

* Minor refactoring of bin/xpimgview.py


*ximpol (0.15.1) - Mon, 25 Jan 2016 16:37:52 -0800*

* astropy.wcs used in ximpol/srcmodel/img.py, and aplpy still used for
  plotting (issue #41).
* Documentation for ximpol/srcmodel/img.py added.


*ximpol (0.15.0) - Mon, 25 Jan 2016 15:57:27 -0800*

* srcmodel config files renamed.
* Point source in the Crab complex sample file dimmer.
* Added option to xpimgview.py to save the image to file.
* Horrible hack in the azimuthal fit to prevent the visibility from going
  negative (issue #34, significantly more work needed).
* Some refactoring and more documentation.
* Radius removed from the xROIModel class, and ROI model for the Crab
  nebula now correctly centered on the right coordinates.


*ximpol (0.14.0) - Fri, 22 Jan 2016 20:54:23 -0800*

* xpobbsim.py generating an output file name based on the source model
  (if not specified).
* Added CMAP mode to xpbin.py


*ximpol (0.13.0) - Fri, 22 Jan 2016 13:58:51 -0800*

* Implemented the infrastructure for multiple source in ROI

  
*ximpol (0.12.1) - Fri, 22 Jan 2016 06:44:01 -0800*

* Bug fix in srcmodel/source.py.


*ximpol (0.12.0) - Thu, 21 Jan 2016 16:35:14 -0800*

* First implementation of extended sources.


*ximpol (0.11.1) - Wed, 20 Jan 2016 16:57:24 -0800*

* Minor addition to the doc.


*ximpol (0.11.0) - Wed, 20 Jan 2016 15:43:39 -0800*

* load_irf moved from bin/xpobssim.py to irf/__init__.py, so that it can be
  reused.
* Unit test for IRF plotting added (issue #30).
* Some documentation for the IRFs added.


*ximpol (0.10.1) - Tue, 19 Jan 2016 16:41:33 -0800*

* More documentation and unit tests.


*ximpol (0.10.0) - Tue, 19 Jan 2016 14:45:50 -0800*

* Added math support in the sphinx config file.
* Major refactoring of the classes related to the modulation factor (issue #28).
* More unit tests added.
* More documentation added.


*ximpol (0.9.1) - Sat, 16 Jan 2016 07:17:52 -0800*

* All unit tests fixed (issue #26).


*ximpol (0.9.0) - Fri, 15 Jan 2016 16:34:58 -0800*

* IRFs extended ("by hand") down below 1 keV (need to do it properly, see
  issue #19).
* A couple of subtle bug fixes in the energy dispersion (see issues #21 and
  #22).
* First version that allows to recover the spectral parameters in XSPEC.


*ximpol (0.8.0) - Fri, 15 Jan 2016 11:53:01 -0800*

* Obsolete files removed, and some name refactoring.
* xpbin.py created.
* All figures from unit tests moved to doc/figures.
* More unit tests.
* Event times in xpobbsim sorted.
* Spectral analysis in xspec added.


*ximpol (0.7.0) - Thu, 14 Jan 2016 15:15:44 -0800*

* Modulation factor generator returning angles in degrees.
* Unit test for the modulation factor classes added.
* Source configuration moved out of xpobsim.py
* Folder srcmodel/config created.
* Added optimization step for the x grid in
  xInterpolatedBivariateSplineLinear.build_vppf() (issue #18).


*ximpol (0.6.3) - Wed, 13 Jan 2016 16:16:38 -0800*

* .travis.yml file tweaked to add display support for matplotlib.


*ximpol (0.6.2) - Wed, 13 Jan 2016 16:11:55 -0800*

* One more unit test added.


*ximpol (0.6.1) - Wed, 13 Jan 2016 15:38:20 -0800*

* Parameter tweak in the xEnergyDispersionMatric class.
* Added unit test for the xCountSpectrum class, with inline images.
* One unit test relaxed.


*ximpol (0.6.0) - Wed, 13 Jan 2016 12:13:06 -0800*

* Number of XIPE energy channels changed from 1024 to 256 and IRFs
  regenerated.
* Removed all the hard-coded values for the number of energy channels
  (issue #13).
* xEnergyDispersionMatrix now inheriting from xUnivariateAuxGenerator (i.e.,
  it has facilities to throw random numbers.)
* Down-sampling mechanism implemented for the xEnergyDispersionMatrix class
  on the energy axis to streamline performance.


*ximpol (0.5.0) - Tue, 12 Jan 2016 15:24:17 -0800*

* A couple of bug fixes in the irf.mrf module.
* Major xpobbsim refactoring using all the new classes.


*ximpol (0.4.2) - Mon, 11 Jan 2016 07:08:21 -0800*

* Minor refactoring.


*ximpol (0.4.1) - Sun, 10 Jan 2016 08:01:03 -0800*

* Grid optimization for the spline definition implemented (issue #15).
* Small application for visualizing an event file (xpevtview.py) created,
  and plotting stuff moved out of xpobbsim.


*ximpol (0.4.0) - Sat, 09 Jan 2016 10:17:52 -0800*

* New module ximpol.core.rand created (issue #16).
* Major rework and speed up of the provisional observation simulator (event
  loop removed).
* New event list classe in.
* Some cleanup.


*ximpol (0.3.1) - Thu, 07 Jan 2016 16:36:04 -0800*

* Added PSF classes, with facility to draw random numbers.


*ximpol (0.3.0) - Thu, 07 Jan 2016 13:53:07 -0800*

* Added make_ppf to the spline base class.
* Some improvement in the plotting facility for the energy dispersion.
* Added unit tests for the irf classes.
* Removed the xmin and xmax arguments from the constructor of all the spline
  classes, since the integral() method does not understand extrapolations and
  having spurious values outside the array ranges was causing troubles.
  (Note the splines can still be extrapolates in the evaluation.)
* Added facilities for normalization, cdf and ppf in the univariate spline
  base class.
* xmerge() method of the base univariate spline class removed in favor of
  numpy.union1d()


*ximpol (0.2.1) - Thu, 07 Jan 2016 06:57:12 -0800*

* First full implementation of the energy dispersion.


*ximpol (0.2.0) - Wed, 06 Jan 2016 15:56:38 -0800*

* Refactoring of the core.spline module, and plotting functionalities added.
* Unit tests for the utils.os_ module added.
* Initial import of the utils.matplotlib_ configuration module.
* Added xEffectiveArea class to irf.arf.
* Added xModulation factor class to mrf.arf.
* bin/xpirfview.py application created (issue #7).


*ximpol (0.1.2) - Tue, 05 Jan 2016 08:34:30 -0800*

* Minor changes.
  

*ximpol (0.1.1) - Tue, 05 Jan 2016 07:05:43 -0800*

* Minor refactoring of the irf specifications, with the OGIP part now included
  in ximpol.irf.base
* Some documentation added to the irf classes.


*ximpol (0.1.0) - Mon, 04 Jan 2016 16:15:30 -0800*

* setup.py file added (issue #11).
* release folder renamed as tools.
* ximpol.__logging__ module moved to ximpol.utils.logging_ (issue #8).
  Note we use the trailing undescore to avoid name conflicts with the
  correponding module from the standard library.)
* ximpol.__utils__ module splitted into ximpol.utils.os_ and
  ximpol.utils.system_ (issue #8).
* Code to create the instrument response functions moved to detector.xipe.
* New spline code used when generating the response functions and old
  xFunction1d classes removed (issue #3).
* fileio folder removed.
* Using the astropy facilities to generate the fits headers (issue #4).


*ximpol (0.0.16) - Sun, 03 Jan 2016 14:31:56 -0800*

* ximpol is now linked to Travis CI, and the build output is shown and linked
  from the main github page.


*ximpol (0.0.15) - Sat, 02 Jan 2016 07:19:39 -0800*

* xChrono class moved to utils.profile. Documentation and unit tests in place.


*ximpol (0.0.14) - Sat, 02 Jan 2016 06:59:19 -0800*

* Minor formatting fix.


*ximpol (0.0.13) - Sat, 02 Jan 2016 06:56:54 -0800*

* Added a makefile for the unit tests, and some more documentation about it.


*ximpol (0.0.12) - Fri, 01 Jan 2016 07:51:56 -0800*

* Some more edits and additions to the documentation.
* Module core.xInterpolatedUnivariateSpline moved to core.spline.
* __package__.py removed, and content moved to ximol.__init__.py, with all
  imports changed accordingly (issue #10).
* Code to be executed in __main__ moved from test() to main() in all modules
  (since the test code will be in the form of unit tests).


*ximpol (0.0.11) - Thu, 31 Dec 2015 17:19:37 -0800*

* Started migrating the documentation from the github wiki to the rst sphinx
  files, and added more stuff.


*ximpol (0.0.10) - Wed, 30 Dec 2015 07:53:08 -0800*

* Bug fix in the release script (hopefully).

  
*ximpol (0.0.9) - Wed, 30 Dec 2015 07:48:26 -0800*

* Major folder restructuring to make the layout compatible with
  `Read the Docs <https://readthedocs.org/>`_.
* Documentation effort started (issue #1).
* Suite of unit tests started (issue #4).
* These release notes moved to a .rst file (issue #12).
* utils.xFunction1d being replaced by core.xInterpolatedUnivariateSpline


*ximpol (0.0.8) - Mon, 28 Dec 2015 06:29:54 -0800*  

* Added script to generate the rmf file. Still not working perfectly.
* Some folder refactoring.


*ximpol (0.0.7) - Fri, 11 Dec 2015 13:33:49 -0800*
  
* Removed the srcmodel/yaml folder and all the associated parser classes.

  
*ximpol (0.0.6) - Fri, 11 Dec 2015 06:39:21 -0800*
  
* Many minor changes.
* First stab at a parser for the source model.
* FITS images of some sources added, along with a small visualization script.
* Added a script that generates the header for the mrf file.
* Added a script to generate the .mrf file based on the ascii table provided.


*ximpol (0.0.5) - Tue, 08 Dec 2015 11:41:24 -0800*
  
* Small fix in the .arf XIPE file.


*ximpol (0.0.4) - Tue, 08 Dec 2015 11:33:40 -0800*
  
* Added a first stab at the effective area table definition.
* Added ascii data files for the XIPE IRFs (as in the proposal).
* Script to generate the .arf file for XIPE based on the ascii table.
* Added a general-purpose one-dimensional function class.


*ximpol (0.0.3) - Fri, 04 Dec 2015 12:11:49 -0800*
  
* Changed thge release note because I was cheating...


*ximpol (0.0.2) - Fri, 04 Dec 2015 12:05:42 -0800*
  
* Folder structure created


*ximpol (0.0.1) - Fri, 04 Dec 2015 06:39:19 -0800*
  
* Initial setup of the github repository.
