Code development
================

Up and running with github
--------------------------

`git <http://git-scm.com/>`_ is a distributed version control system and
`github <https://github.com/>`_ is the web hosting service that we use to
develop and maintain `ximpol <https://github.com/lucabaldini/ximpol>`_.

You will need to register to github and ask
`Luca Baldini <mailto:luca.baldini@pi.infn.it>`_ to be granted write-access to
the repository. `Here <http://git-scm.com/doc>`_ is the entry point for the git
documentation, in case you want to have a feeling of what git is doing and how
you use it.

Mind that, in order to be able to push back changes to the remote repository
you will need to tell git on your machine who you are, i.e.:

>>> git config --global user.name "Your Name"
>>> git config --global user.email you@example.com


Cloning the repository - GNU/Linux and Mac OS
---------------------------------------------

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



Coding conventions
------------------


Documenting code
----------------

We use `sphinx <http://sphinx-doc.org/#>`_ to generate the ximpol
`documentation <http://ximpol.readthedocs.org/en/latest/index.html>`_.
