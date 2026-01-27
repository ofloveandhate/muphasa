🛠 Installation
****************************

Installing `muphasa`
-------------------------

Muphasa can be installed with `pip`.  The build system is `setuptools`, and this library requires `cython`.  

You need the following things in order to install muphasa:

* numpy
* cython
* a c++ compiler, that matches the version of python installed
* boost, matching the c++ compiler

Getting the environment right is easiest using a package manager.  We have found `micromamba` a successful tool.

1. Install `micromamba`.  Get it activated in your shell, and make an environment.  Activate that environment.
2. Install the dependencies.  `mamba install numpy cython python boost`.  Depending on other things you have installed, you may need to include `g++`.
3. Install muphasa.  Move to the directory and `pip install .`.

If you wish to use tools like `ipython`, `jupyter`, etc, install those through the package manager to ensure that all the various things line up.


Building the documentation
-----------------------------

The documentation is built using `sphinx`.  Install it via the package manager, `mamba install sphinx sphinxcontrib-bibtex`

Either you can fully install `mph` to your environment, then build the documentation from the installed version, or you can build from the checked out source:

Build from checked out source:

1. `pip install -e .`  This builds the compiled cython part of muphasa locally.
2. from `docs/`, run `make html`.  This will use the local version.
