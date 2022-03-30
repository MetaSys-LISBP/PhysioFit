Frequently Asked Questions (FAQ)
================================

What units should be provided?
------------------------------

Input data (biomass concentration, metabolites concentrations, and time) can be provided to PhysioFit 2.0 using any unit. To ensure numerical stability, we recommand to use units as close to unity as
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

What is a Khi2 test?
------------------------------------------------------------------

XXX

My data hasn't been correctly fitted. Why?
------------------------------------------------------------------

There are many reasons why a fit can be incorrect. In most cases, it is because the run parameters have to be tweaked to
help PhysioFit get on the right track. For more info on the run parameters and how they affect the fitting process,
please refer to section :ref:`physiofit parameters`.
Importantly, a key asumption of the flux calculation method is that cells are in metabolic (pseudo) steady-state, i.e. with constant fluxes during the whole experiment. If this is not the case (e.g. cells are continuously adapting to their environment and fluxes change over time), PhysioFit will not be able to fit the data satisfactorily.

