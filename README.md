# PhysioFit

[![PyPI version](https://badge.fury.io/py/physiofit.svg)](https://badge.fury.io/py/physiofit)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/physiofit.svg)](https://pypi.python.org/pypi/physiofit/)
[![Documentation Status](https://readthedocs.org/projects/physiofit/badge/?version=latest)](http://physiofit.readthedocs.io/?badge=latest)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/physiofit/README.html)


## What is PhysioFit?
**PhysioFit is a scientific tool designed to quantify cell growth parameters and uptake & production fluxes**

Fluxes are estimated using mathematical models by fitting time-course measurements of the concentration of
cells and extracellular substrates and products. PhysioFit is shipped with some common growth models, and
additional tailor-made models can be implemented by users.

It is one of the routine tools that we use at the [MetaSys team](https://www.toulouse-biotechnology-institute.fr/en/poles/equipe-metasys/) 
and [MetaToul platform](https://www.metabohub.fr/home.html) in functional studies of metabolic systems.

The code is open-source, and available under a GPLv3 license. Additional information can be found in the following 
[publication](https://doi.org/10.1101/2023.10.12.561695).

Detailed documentation can be found online at Read the Docs 
([https://physiofit.readthedocs.io/](https://physiofit.readthedocs.io/)).

## Key features

   * **calculation of growth rate and extracellular (uptake and production) fluxes**,
   * **a set of steady-state and dynamic models** are shipped with PhysioFit,
   * **tailor-made models** can be constructed by users,
   * Monte-Carlo sensitivity analysis to **estimate the precision on the calculated fluxes**,
   * **evaluation of the goodness of fit and visual inspection of the fitted curves**,
   * shipped as a **library** with both a **graphical** and a **command line** interface,
   * **open-source, free and easy to install** everywhere where Python 3 and pip run,
   * **biologist-friendly**.

## Quick-start
PhysioFit requires Python 3.9 or higher and run on all platforms.
Please check [the documentation](https://physiofit.readthedocs.io/en/latest/quickstart.html) for complete
installation and usage instructions.

Use `pip` to **install PhysioFit from PyPi**:

```bash
$ pip install physiofit
```

Then, start the graphical interface with:

```bash
$ physiofit
```

PhysioFit is also available directly from command-line and as a Python library.

## Bug and feature requests
If you have an idea on how we could improve PhysioFit please submit a new *issue*
to [our GitHub issue tracker](https://github.com/MetaSys-LISBP/PhysioFit/issues).


## Developers guide
### Contributions
Contributions are very welcome! :heart:

Please work on your own fork,
follow [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide,
and make sure you pass all the tests before a pull request.

### Local install with pip
In development mode, do a `pip install -e /path/to/PhysioFit` to install
locally the development version.

### Build the documentation locally
Build the HTML documentation with:

```bash
$ cd doc
$ make html
```

The PDF documentation can be built locally by replacing `html` by `latexpdf`
in the command above. You will need a recent latex installation.

## How to cite
PhysioFit: a software to quantify cell growth parameters and extracellular fluxes.
Le Grégam L., Guitton Y., Bellvert F., Jourdan F., Portais J.C., Millard P.
bioRxiv preprint, [doi: 10.1101/2023.10.12.561695](https://doi.org/10.1101/2023.10.12.561695)

## Authors
Loïc Le Grégam, Pierre Millard

## Contact
:email: legregam@insa-toulouse.fr, millard@insa-toulouse.fr
