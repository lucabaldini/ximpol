Installation
============

Prerequisites
-------------

The package is based on the `Python <https://www.python.org/>`_ scripting
language and the `SciPy <http://www.scipy.org/>`_ Python-based ecosystem.
You will need a working Python installation including several different
side packages, most notably:

* `NumPy <http://www.numpy.org/>`_: the fundamental package for scientific
  computing with Python. 
* `SciPy <http://www.scipy.org/>`_: a Python-based ecosystem of open-source
  scientific software. 
* `matplotlib <http://matplotlib.org/>`_: a Python 2D plotting library.
* `Astropy <http://www.astropy.org/>`_: a Python astronomy package (including
  tools to manipulate FITS files).
* `APLpy <https://aplpy.github.io/>`_ (optional): a Python module for plotting
  astronomical imaging data in FITS format.
* `pyregion <http://pyregion.readthedocs.io/en/latest/>`: a Python module to
  parse ds9 region files
* `PyXspec <https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/>`_
  (optional): the Python binding for XSPEC.

Loosely speaking, you should be able to open the Python terminal and execute
the following ``import`` statements with no errors.

>>> import numpy
>>> import scipy
>>> import matplotlib
>>> import astropy
>>> import aplpy
>>> import xspec
>>> import pyregion

(The last two are not strictly necessary to run a simulation, but in the
long run you probably want to have them.) If any of the required packages
fails to import, take a deep breath and fix the issue before you move on.

Using anaconda
--------------

Since the all the python third-party packages that we use are under active
development (and for some of them there have been some non-trivial
interface changes across the most recent version), you do care about the
package versions that you have installed. For instance, we do need the
ability to evaluate spline on multidimensional arrays that was added in
SciPy 0.16 but was not there in SciPy 0.14. The bottomline is that, depending
on the exact version of the packages that you have installed on your system,
ximpol might or might not work correctly.

One possible OS-independent way to get a fully fledged ecosystem in which
ximpol can work is to use `anaconda <https://www.continuum.io/downloads>`_.
You should be able to get up and running, in terms of the pre-requisites to
run ximpol, in a matter of minutes following the instructions in there.

Note that the default anaconda installer does not come with all the packages
that you need (most notably, aplpy is missing), but you can easily install
whetever you need after the fact, e.g.,

.. code-block:: bash
                
    pip install aplpy
    pip install pyregion

If you are uncertain on what to do, you might want to try anaconda first,
as this might be the less painful solution.

For XSPEC and the corresponding Python bindings, refer to the
`HEASOFT download page <http://heasarc.nasa.gov/lheasoft/download.html>`_.


More OS-specific resources
--------------------------

There are obviously other ways to install the prerequisites, e.g., by
doing everything manually or use a package manager that your OS makes available.

.. warning:: We need to add specific information for different platforms
             (e.g., Windows and Mac) and GNU/Linux distros. Any help from
             anybody is appreciated!

For GNU/Linux in the Fedora flavor, for instance, you would do something like

.. code-block:: bash

    dnf install numpy scipy python-matplotlib python-astropy APLpy


Downloading the code
--------------------

The easiest (though not the best) way to get the code is by directly
downloading the zip or tar file for the latest tag (or the tag you are
interested in) from the `github release page
<https://github.com/lucabaldini/ximpol/releases>`_. Unzip the archive in
your favorite folder and setup the environment as detailed in the next
section.

If you plan on actively contributing to the software development (as opposed
to just using it) you will need to clone the github repository, as explained
in the page about the code development.

.. warning:: At some point we should provide a working cross-platform
             installation mechanism through the setup.py file, and possibly
             packages/installers for the operating systems we use. Again,
             any help would be appreciated.


Basic environment
-----------------

The only thing you have to do is to make sure that the root folder of the
repository is included in the ``$PYTHONPATH`` environmental variable.
You might also want to add `ximpol/bin` to the ``$PATH`` environmental variable,
so that you have the executables off hand. Here is an example for users of the
Bourne shell (sh, ash, ksh, and bash): 

.. code-block:: bash

    export PYTHONPATH=/data/work/xipe/ximpol:$PYTHONPATH
    export PATH=/data/work/xipe/ximpol/ximpol/bin:$PATH

Loosely speaking, if you can open a Python prompt and do

>>> import ximpol

without getting back an error message like this

>>> import ximpol
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: No module named ximpol
>>> 

again, you should be all set. If not, I am sorry to say, you really do have to
fix this before moving on.
