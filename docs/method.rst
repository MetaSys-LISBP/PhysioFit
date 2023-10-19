Method
===============

General principle
*****************

Fluxes (exchange fluxes - qM\ :sub:`i` - and growth rate - µ -), initial concentrations of species (i.e. biomass - X - and 
metabolites - M_i - at t=0) and possibly other growth parameters (e.g. lag time) are estimated by fitting time-course measurements of
metabolite and biomass concentrations, as detailed in this section.

Flux values provided by PhysioFit correspond the best fit. A global sensitivity analysis (Monte-Carlo approach) is
available to evaluate the precision of the estimated fluxes (mean, median, standard deviation, 95% confidence
intervals), plots are generated for visual inspection of the fitting quality, and a χ² test is performed to assess the
statistical goodness of fit.

.. _method_models:

Models
******

Models are at the heart of the flux calculation approach implemented in PhysioFit. A flux  model contains equations that describe the dynamics of biomass and 
metabolite concentrations (thus allowing simulation of metabolite dynamics) and the list of all parameters (including fluxes) with their 
initial values and bounds (used for optimization and flux calculation). 

Different models are shipped with PhysioFit, and tailor-made models can be provided by users, as detailed in the :doc:`models` section.

.. _optimization_process:

Flux calculation process
************************

A model containing all metabolites provided in the input data is constructed by PhysioFit. 
Model parameters (such as fluxes, growth rate, and initial concentrations of biomass and metabolites) are estimated by fitting experimental metabolite and biomass dynamics, where 
the following cost function is minimized:

.. image:: _static/equations/eq10.png

where *sim* is the simulated data, *meas* denotes measurements, and *weight* is a weighting factor (typically, one
standard deviation on measurements).

For this optimization step, PhysioFit uses the Scipy's Differential evolution method to approximate the solution, and the best solution is polished using the L-BFGS-B method (see
`scipy.optimize <https://docs.scipy.org/doc/scipy/reference/optimize.html>`_ for more information on the optimization
method).

Sensitivity analysis
*********************

To quantify the standard deviation and 95% confidence intervals on the fits of
experimental data & on the estimated parameters (including fluxes), we use a Monte Carlo analysis. Briefly, PhysioFit generates a
number of noisy datasets from the experimental data (i.e option "number of iterations" in the GUI) and calculate fluxes and other growth 
parameters for each of the synthetic datasets. This enables PhysioFit to perform statistical analysis and estimate the standard deviation and confidence interval for 
each parameter. We recommend to always run flux calculations with sensitivity analyses.

