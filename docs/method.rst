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


.. _sensitivity_analysis:

Sensitivity analysis
*********************

To determine the precision on the fit and on the estimated parameters,
PhysioFit performs a Monte Carlo analysis. Briefly,
PhysioFit generates several
datasets by adding noise to the dynamics simulated from the best fit, and
calculates fluxes and other growth
parameters for each of these synthetic datasets. This enables PhysioFit to
compute statistics (mean, median, standard deviation and 95 % confidence
interval) for each parameter. We recommend always running
a sensitivity analysis when using PhysioFit.


..  _`chi2 test`:

Goodness-of-fit evaluation
**************************

PhysioFit performs a χ² test to assess the goodness of fit. A χ² test
describes how well a model fits a set of observations. Measures of
goodness of fit typically summarize the discrepancy between observed values
and the values expected under the model used in PhysioFit (see
:ref:`optimization_process`). It is calculated as the sum of differences
between measured and simulated values, each squared and divided by the
simulated value.
A good fit corresponds to small differences between measured and simulated
values, thereby the χ² value is low. In contrast, a bad fit corresponds to
large differences between simulations and measurements, and the χ² value is
high.

The resulting χ² value can then be compared with a χ² distribution to
determine the goodness of fit. The p-value of one-tail χ² test is calculated
by PhysioFit from the best fit and is given in the log file (have a look to
the :doc:`usage` section). A p-value close to 0 means poor fitting, and a
p-value close to 1 means good fitting (keeping in mind that a p-value very
close to 1 can be an evidence that standard deviations might be
overestimated). A p-value between 0.95 and 1 means the model fits the data
good enough with respect to the standard deviations provided (at a 95%
confidence level). PhysioFit provides an explicit message stating whether
the flux data are satisfactorily fitted or not (at a 95% confidence interval).

Model comparison
***********************************

PhysioFit calculates the Aikake Information Criterion (AIC)
to help users select the most appropriate model for their data.
The AIC is a measure of the explanatory power of a model with respect to a
given set of data. The model with the lowest AIC value is
thus considered the best model.

The AIC is calculated as follows:

.. math::

    AIC = 2k + n \ln(\frac{RSS}{n})

where :math:`k` is the number of parameters in the model (plus 1), :math:`n` is the
number of data points, and :math:`RSS` is the residual sum of squares. For
small sample sizes (typically < 40 data points>), it is recommended to use 
the AICc (corrected AIC), which is
calculated as follows:

.. math::

    AICc = AIC + \frac{2k(k+1)}{n-k-1}

In practice, because the AICc approximates the AIC for large sample sizes,
it's often advised that AICc be used as default.

The AIC
value can be found in the statistical output file or directly in
the graphical user interface. 

To identify the best model, different candidate models that differ in terms
of structure or complexity can be used
to fit the data and then compared based on their AIC. The model
with the lowest AIC value is considered the best-fitting model among
the candidates and should thus used to fit the dataset. However, it is crucial to consider the differences
in AIC values between models, as models with low ΔAIC values (typically < 2)
are considered to have similar support from the data. 

More information on the AIC can be found in the `original paper by Akaike
(1974) <https://gwern.net/doc/statistics/decision/1998-akaike.pdf>`_ or `in
this practical guide <https://doi.org/10
.1007/s00265-010-1037-6>`_.

