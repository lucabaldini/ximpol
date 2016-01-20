Installation
============

Prerequisites
-------------

The package is based on the `Python <https://www.python.org/>`_ scripting
language and the `SciPy <http://www.scipy.org/>`_ Python-based ecosystem.
You will need a working Python installation including the SciPy library and
`matplotlib <http://matplotlib.org/>`_.

ximpol also uses `Astropy <http://www.astropy.org/>`_ for a few things (most
notably the
`astropy.io.fits <http://docs.astropy.org/en/stable/io/fits/index.html>`_
module for manipulating FITS files), so you'll need that as well.


Cloning the repository
----------------------

(This will work on GNU/Linux and MAC. We should provide instructions for
Windows as well.)

Cloning the reprository is as easy as typing

>>> git clone git@github.com:lucabaldini/ximpol.git

(mind this will create a ximpol folder in the current directory, so cd to the
appropriate place in your filesystem first). If you get an error message along
the lines of

>>> Permission denied (publickey).
>>> fatal: Could not read from remote repository.
>>> Please make sure you have the correct access rights and the repository exists.

that simply means that you have to exchange you public SSH key with the server.
In order to do that, click on your github profile icon on the top-right of the
github webpage, select "Edit profile", "SSH keys" (on the left-hand side
"Personal setting" menu), "Add SSH key" and paste in the form the content of
the local (i.e. on the machine you are cloning the repository into)
`~/.ssh/id_rsa.pub` file.

If you don't have a public ssh key, you can generate it by typing

>>> ssh-keygen

(press ENTER a couple of times and here is you public key in
`~/.ssh/id_rsa.pub`)


Basic environment
-----------------

You have to make sure that the root folder of the repository is included in
the `$PYTHONPATH` environmental variable.

You might also want to add `ximpol/bin` to the `$PATH` environmental variable,
so that you have the executables off hand. Here is an example:

>>> export PYTHONPATH=/data/work/xipe/ximpol:$PYTHONPATH
>>> export PATH=/data/work/xipe/ximpol/ximpol/bin:$PATH

