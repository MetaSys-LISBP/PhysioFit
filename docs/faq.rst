Frequently asked questions (FAQ)
================================

How are fluxes calculated?
------------------------------------------------------------------

We provide details on the flux calculation approach implemented in PhysioFit in the :doc:`method` section.

What units should be provided?
------------------------------

Input data (biomass concentration, metabolites concentrations, and time) can be provided to PhysioFit using any unit. To ensure numerical stability, we recommand to use units for which values are as close to unity as
possible (e.g. 3 mM instead of 3e3 µM). Importantly, units of the extimated fluxes depend on units of time and metabolites and biomass concentrations 
(see below).

What are the flux units?
------------------------

Flux units depend on the units of time and concentrations (of biomass and metabolites). For more info, please refer to the
:doc:`usage` section.

An error has been raised. What should I do?
-------------------------------------------

The first thing to do is to read the error message which might contain information on how to resolve it. If not, check the FAQ
section (yes, this one) to see if the error has been explained in more depth. If the error persists or if you do not
understand the error, please post it in the "issues" section on `GitHub
<https://github.com/MetaSys-LISBP/PhysioFit/issues>`_. We will try to respond as quickly as possible to solve your problem.

What parameter values should I use?
------------------------------------------------------------------

TO BE IMPLEMENTED

How can I check if my data have been fitted correctly?
------------------------------------------------------------------

TO BE IMPLEMENTED

What is a χ² test?
------------------------------------------------------------------

Pearson's χ² test is a measure of goodness of fit. It is calculated as the sum of differences between measured and
simulated values, each squared and divided by the simulated value. The resulting value can then be compared with a khi2
distribution to determine the goodness of fit. A good fit corresponds to small differences between measured and
simulated values, thereby the χ² value is low. In contrast, a bad fit corresponds to large differences, and the χ² value
is high. Once fluxes are determined, PhysioFit calculates the p-value of one-tail χ² test, which is given in the log
file (have a look to the :doc:`usage` section). Value close to 0 (e.g. under 0.05) means poor fitting, and value between
0.95 and 1 means good fitting (according to the provided weights, which correspond typically to the standard deviations
on measurements). Value very close to 1 can be an evidence that weights are overestimated.

My data hasn't been correctly fitted. Why?
------------------------------------------------------------------

There are many reasons why a fit can be incorrect. In most cases, it is because the run parameters have to be tweaked to
help PhysioFit get on the right track. For more info on the run parameters and how they affect the fitting process,
please refer to section :ref:`physiofit parameters`.
Importantly, a key asumption of the flux calculation method is that cells are in metabolic (pseudo) steady-state, i.e. with constant fluxes during the whole experiment. If this is not the case (e.g. cells are continuously adapting to their environment and fluxes change over time), PhysioFit will not be able to fit the data satisfactorily.

I would like a new feature.
------------------------------------------------------------------

We would be glad to improve PhysioFit. Please get in touch with us <https://github.com/MetaSys-LISBP/PhysioFit/issues>`_ so we could discuss your problem.
