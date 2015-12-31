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



Coding guidelines
-----------------

Though we'll never be able to follow any set of coding conventions religiously,
`PEP 0008 <https://www.python.org/dev/peps/pep-0008/>`_ is our starting point.
A short recap of the most salient guidelines:

* Use 4 spaces for indentation level (no TABS).
* Limit all lines to 79 characters.
* Surround top-level function and class definitions with two blank lines.
  Method definitions inside a class are surrounded by a single blank line.
  Use blank lines in functions, sparingly, to indicate logical sections.
* Use one import per line, right at the top of the module.
* Use single quote characters for strings and double quotes characters for 
  triple-quoted strings.
* Avoid extraneous white spaces, and especially avoid more than one space
  around an assigment.
* Don't use spaces around the `=` sign when used to indicate a keyword argument
  or a default parameter value.
* Modules should have short, all-lowercase names.
* Class names should normally use the CapWords convention (for ximpol starting
  with a `x`).
* Function and member names should be lowercase, with words separated by
  underscores as necessary to improve readability.
* Constants are usually defined on a module level and written in all capital
  letters with underscores separating words.
* Always use a `def` statement instead of an assignment statement that binds a
  `lambda` expression directly to an identifier. 

An example module, illustrating the basic guidelines, is available
`here <https://github.com/lucabaldini/ximpol/tree/master/doc/codestyle.py>`_.


Documenting code
----------------

We use `sphinx <http://sphinx-doc.org/#>`_ to generate the ximpol
`documentation <http://ximpol.readthedocs.org/en/latest/index.html>`_.
