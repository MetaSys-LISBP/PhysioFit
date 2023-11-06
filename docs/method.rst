Method
===============

Overview
*********

Fluxes (exchange fluxes of a metabolite M\ :sub:`i`, qM\ :sub:`i` ; growth rate, µ), initial concentrations of species (biomass, X ; 
metabolites, M\ :sub:`i`) and possibly other growth parameters (e.g. lag time) are estimated by fitting time-course measurements of
metabolite and biomass concentrations, as detailed below.

Flux values provided by PhysioFit correspond the best fit. A global sensitivity analysis (Monte-Carlo approach) is
available to evaluate the precision of the estimated fluxes (mean, median, standard deviation, 95% confidence
intervals), plots are generated for visual inspection of the fitting quality, and a χ² test is performed to assess the
statistical goodness of fit.

.. _method_models:

Models
******

Models are at the heart of the flux calculation approach implemented in PhysioFit. A flux model contains i) equations that describe the dynamics of biomass and 
metabolite concentrations as function of different parameters (used to simulate time-course metabolite concentrations) and ii) the list of all parameters (including fluxes) with their 
(default) initial values and bounds (used for flux calculation). 

Different models are shipped with PhysioFit, and tailor-made models can be provided by users, as detailed in the :doc:`models` section.

.. _optimization_process:

Flux calculation
************************

First, PhysioFit construct a model that used to simulate the dynamics of the concentration of biomass and metabolites (substrates and products) provided in the input data. 
Model parameters (such as fluxes, growth rate, and initial concentrations of biomass and metabolites) are then estimated by fitting experimental metabolite and biomass dynamics. PhysioFit 
minimizes the following cost function:

.. math:: residuum = \sum_{i} (\dfrac{sim_{i}-meas_{i}}{sd_{i}})^2

where :math:`sim` is the simulated data, :math:`meas` denotes measurements, and :math:`sd` is the 
standard deviation on measurements.

For this optimization step, PhysioFit uses the Scipy's Differential evolution method to approximate the solution, 
and the best solution is polished using the L-BFGS-B method (see
`scipy.optimize <https://docs.scipy.org/doc/scipy/reference/optimize.html>`_ for more information on the optimization
process).

Goodness-of-fit evaluation
**************************

PhysioFit performs a χ² test to assess the goodness of fit. Have a look at the :doc:`faq` section for 
more details on the interpretation of the khi2 test results.

Sensitivity analysis
*********************

To determine the precision on the fit and on the estimated parameters (including fluxes), PhysioFit performs a Monte Carlo analysis. Briefly, PhysioFit generates several 
datasets by adding noise to the dynamics simulated from the best fit, and calculated fluxes and other growth 
parameters for each of these synthetic datasets. This enables PhysioFit to compute statistics (mean, median, standard deviation and 95% confidence interval) for 
each parameter (including fluxes). We recommend always running a sensitivity analysis when using PhysioFit.

