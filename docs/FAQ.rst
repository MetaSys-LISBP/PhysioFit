Frequently Asked Questions (FAQ)
================================

What units should be provided?
------------------------------

Any units can be provided to PhysioFit 2.0. It is however recommended to use the unit that is as close to unity as
possible (e.g. 3 mM instead of 3e3 µM). Of note, units of the extimated fluxes depend on units of time and metabolites and biomass concentrations 
(see below).

What are the flux units?
------------------------

Flux units depend on the units of time and concentrations. For more info, please refer to the **Warning** block in the
:doc:`usage` section.

An error has been raised. What should I do?
-------------------------------------------

The first thing to do is to read the error message which might contain info on how to resolve it. If not, check the FAQ
section (yes, this one) to see if the error has been explained in more depth. If the error persists or if you do not
understand the error, please  post it in the "issues" section on `GitHub
<https://github.com/MetaSys-LISBP/PhysioFit/issues>`_. We will try to respond as quickly as possible with a fix.

After flux calculation, my data hasn't been correctly fitted. Why?
------------------------------------------------------------------

There are many reasons why a fit can be incorrect. In most cases, it is because the run parameters have to be tweaked to
help PhysioFit get on the right track. For more info on the run parameters and how they affect the fitting process,
please refer to the :ref:`physiofit parameters`.
