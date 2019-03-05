Repository for QuTiP documentation
==================================

This repository contains the source files for the QuTiP documentation.

For pre-built documentation, see http://www.qutip.org/documentation.html

Build requirements
------------------

* Sphinx: http://sphinx-doc.org/
* sphinx_rtd_theme: https://sphinx-rtd-theme.readthedocs.io/
* numpydoc: https://numpydoc.readthedocs.io/en/latest/
* sphinxcontrib-bibtex: https://github.com/mcmtroffaes/sphinxcontrib-bibtex
* LaTeX and pdflatex

Install using pip:
    
    $ pip install sphinx numpydoc sphinx_rtd_theme sphinxcontrib-bibtex

2019-03-05: Successful building using:

* sphinx v1.8.4
* sphinx_rtd_theme v0.4.3
* numpydoc v0.8.0
* 

Build
-----
To build the documentation on Linux or OS X run:

    $ make html
