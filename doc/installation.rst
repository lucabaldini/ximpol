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

(The last two are not strictly necessary to run a simulation, but in the
long run you probably want to have them.) If any of the required packages
fails to import, take a deep breath and fix the issue before you move on.

.. warning:: We need to add specific information for different platforms
             (e.g., Windows and Mac) and GNU/Linux distros. Any help from
             anybody is appreciated!

How you actually go about making sure that all the packages are correctly
installed depends on the operating system you're using. For GNU/Linux
in the Fedora flavor, for instance, you would do something like

.. code-block:: bash

    yum install numpy scipy python-matplotlib python-astropy APLpy

For XSPEC and the corresponding Python bindings, refer to the
`HEASOFT download page <http://heasarc.nasa.gov/lheasoft/download.html>`_.


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
