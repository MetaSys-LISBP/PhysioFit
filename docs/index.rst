.. PhysioFit documentation master file, created by
   sphinx-quickstart on Mon Mar  7 14:29:46 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/Logo.png
   :align: left
   :scale: 36%

Welcome to PhysioFit documentation!
=====================================

**PhysioFit is a scientific tool designed to i) quantify exchange (production and consumption) fluxes and ii) cell growth
rate during (batch) cultivations of microorganisms.**

Fluxes are estimated from time-course measurements of extracellular metabolites and biomass concentrations. PhysioFit has been designed to 
calculate fluxes in batch experiments, assuming cells are in metabolic (pseudo) steady-state (i.e. fluxes are constant during the experiment).

**PhysioFit includes the following features:**

   * **calculation of growth rate and extracellular (uptake and production) fluxes**.
   * if cell growth has some **lag** (e.g. due to adaptation to a novel environment), lag can be taken into account and lag time estimated.
   * **non-enzymatic degradation** of some metabolites (e.g. DHA or glutamine) can be estimated and taken into account when calculating fluxes.
   * sensitivity analyses are performed to **estimate the precision of the calculated fluxes**.
   * **evaluation of the goodness of fit and visual inspection of the fitted curves**.
   * shipped as a **library** with both a **graphical** and **command line** interface,
   * open-source, free and easy to install everywhere where Python 3 and pip run,
   * biologist-friendly.

It is one of the routine tools that we use at the `MetaSys team <http://www.toulouse-biotechnology-institute.fr/en/research/molecular-physiology-and-metabolism/metasys.html>`_ and `MetaToul platform <http://www.metatoul.fr>`_ to calculate fluxes.

The code is open-source, and available on `GitHub <https://github.com/MetaSys-LISBP/PhysioFit/>`_ under a :ref:`GPLv3 license <license>`.

This documentation is available on Read the Docs (`https://physiofit.readthedocs.io <https://physiofit.readthedocs.io/>`_)
and can be downloaded as a `PDF file <https://readthedocs.org/projects/physiofit/downloads/pdf/latest/>`_.


.. toctree::
   :maxdepth: 2
   :caption: Usage

   installation.rst
   quickstart.rst
   usage.rst
   method.rst
   cite.rst


.. toctree::
   :maxdepth: 2
   :caption: Miscellaneous

   librairy.rst
   faq.rst
   license.rst

* :ref:`search`
