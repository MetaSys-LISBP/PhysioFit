Frequently asked questions (FAQ)
================================

How are fluxes calculated?
------------------------------------------------------------------

We provide details on the flux calculation approach implemented in PhysioFit in the :doc:`method` section.

What units should be used for input data?
------------------------------

Input data (biomass concentration, metabolites concentrations, and time) can be provided to PhysioFit using any unit. Still, we recommand to use units for which values are as close to unity as
possible to ensure numerical stability (e.g. 3 mM instead of 3e3 µM). Importantly, units of the extimated fluxes depend on units of time and metabolites and biomass concentrations, as detailed below. Additional 
information on units can be found in the :doc:`usage` section.

What are the flux units?
------------------------

Flux units depend on the units of time and concentrations (of biomass and metabolites). For more details, please refer to the
:doc:`usage` section.

An error has been raised. What should I do?
-------------------------------------------

The first thing to do is to read the error message which might contain information on how to resolve it. If not, check the FAQ
section (yes, this one) to see if the error has been explained in more depth. If the error persists or if you do not
understand the error, please post it in the "issues" section on `GitHub
<https://github.com/MetaSys-LISBP/PhysioFit/issues>`_. We will try to respond as quickly as possible to solve your problem.

What parameters values should I use?
------------------------------------------------------------------

Details on PhysioFit parameters can be found in the :doc:`usage` section.

How can I check if my data have been fitted correctly?
------------------------------------------------------------------

The quality of the fit can be evaluated based on:

    * the plots of experimental vs simulated data for the best fit.
    * the χ² statistical test results given in the log file (see below for help on interpreting the results).

What is a χ² test?
------------------------------------------------------------------

Pearson's χ² test is a measure of goodness of fit. It is calculated as the sum of differences between measured and simulated values, each squared and divided by the simulated value. 
A good fit corresponds to small differences between measured and simulated values, thereby the χ² value is low. In contrast, a bad fit corresponds to large differences, and the χ² value is high. 

The resulting χ² value can then be compared with a χ² distribution to determine the goodness of fit. The p-value of one-tail χ² test is calculated by PhysioFit from the best fit and is given in the log file (have a look to the :doc:`usage` section). Value close to 0 means poor fitting, and value close to 1 means good fitting (keeping in mind that a p-value very close to 1 can be an evidence that weights might be overestimated). A 
p-value should between 0.95 and 1 means the model fits the data good enough with respect to the standard deviations provided as weights (at level of 95% confidence).

My data hasn't been correctly fitted. Why?
------------------------------------------------------------------

The main reason to explain a bad fit is that a key asumption of the flux calculation method is that cells are in metabolic (pseudo) steady-state, i.e. with constant fluxes during the whole experiment. If this is not the case (e.g. cells are continuously adapting to their environment and fluxes change over time), PhysioFit will not be able to fit the data satisfactorily.

In rare situations, it may be because the run parameters have to be tweaked to
help PhysioFit get on the right track, which results in obviously aberrant fits (e.g. with flat time-course profiles for all metabolites). For more info on the run parameters and how they affect the fitting process,
please refer to section :ref:`physiofit parameters`.

I would like a new feature.
------------------------------------------------------------------

We would be glad to improve PhysioFit. Please get in touch with us <https://github.com/MetaSys-LISBP/PhysioFit/issues>`_ so we could discuss your problem.
