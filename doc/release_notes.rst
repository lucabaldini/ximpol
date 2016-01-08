Release notes
=============

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
