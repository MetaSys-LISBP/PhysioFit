# PhysioFit

[![PyPI version](https://badge.fury.io/py/physiofit.svg)](https://badge.fury.io/py/physiofit)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/physiofit.svg)](https://pypi.python.org/pypi/physiofit/)
[![Documentation Status](https://readthedocs.org/projects/physiofit/badge/?version=latest)](http://physiofit.readthedocs.io/?badge=latest)


## What is PhysioFit?
PhysioFit is a scientific tool designed to i) quantify exchange (production and consumption) fluxes and  ii) cell growth rate during (batch) cultivations of microorganisms. Fluxes are estimated from time-course measurements of extracellular metabolites and biomass concentrations. An important assumption is that cells are in metabolic (pseudo) steady-state.

It is one of the routine tools that we use at the [MetaSys team](http://www.lisbp.fr/en/research/molecular-physiology-and-metabolism/metasys.html) and [MetaToul platform](http://www.metatoul.fr) in functional studies of metabolic systems.

The code is open-source, and available under a GPLv3 license. Additional information can be found in the following [publication](https://doi.org/10.1128/aem.00768-19).

Detailed documentation can be found online at Read the Docs ([https://physiofit.readthedocs.io/](https://physiofit.readthedocs.io/)).

## Key features
* **calculation of growth rate and extracellular (uptake and production) fluxes**.
* **lag** before growth (e.g. due to adaptation to a novel environment) can be taken into account and estimated.
* **non-enzymatic degradation** of some metabolites (e.g. DHA or glutamine) can be estimated and taken into account when calculating exchange fluxes.
* **sensitivity analyses are performed to estimate the precision of the calculated fluxes**.
* **evaluation of the goodness of fit and visual inspection of the fitted curves**.
* shipped as a **library** with both a **graphical** and **command line** interface,
* open-source, free and easy to install everywhere where Python 3 and pip run,
* biologist-friendly.

## Quick-start
PhysioFit requires Python 3.6 or higher and run on all plate-forms.
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
Peiro C., Millard P., de Simone A., Cahoreau E., Peyriga L., Enjalbert B., and Heux S. Chemical and metabolic controls on dihydroxyacetone metabolism lead to suboptimal growth of *Escherichia coli*. Appl Environ Microbiol, 2019, [doi: 10.1128/AEM.00768-19](https://doi.org/10.1128/AEM.00768-19)

## Authors
Loïc Le Grégam, Pierre Millard

## Contact
:email: legregam@insa-toulouse.fr, millard@insa-toulouse.fr
