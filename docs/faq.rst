Frequently asked questions (FAQ)
================================

How are fluxes calculated?
------------------------------------------------------------------

We provide details on the flux calculation approach implemented in PhysioFit in the :doc:`models` section.

How many measurements should I use to calculate fluxes?
------------------------------------------------------------------

As in any model-based fitting procedure, more data means more accurate and precise flux estimates. We recommend using
at least 6 time points, which should provide reliable and meaningful estimates in most situations.

Still, the exact answer to this question strongly depends on the uptake/production/growth rates of your (micro)organism
in the conditions you are investigating, on the sampling time interval, on the questions you are addressing, and on
many other parameters! You can make some tests by calculating fluxes from (published or theoretical) datasets similar
to those you have in mind.

Can I calculate fluxes in case of missing values?
------------------------------------------------------------------

Yes, fluxes can still be calculated if some measurement(s) are missing. In this case, let empty the corresponding field
of the input data file.

..  _`conc units`:

What units should be used for input data?
-----------------------------------------

Input data (biomass concentration, metabolites concentrations, and time) can be provided to PhysioFit using any unit.
Still, we recommand to use units for which values are as close to unity as possible to ensure numerical stability (e.g.
3 mM instead of 3e3 µM). Importantly, units of the estimated fluxes depend on units of time and metabolites and biomass
concentrations. The concentration of different metabolites can be provided using different units, but a single unit
must be used for all measurements of a given metabolite.

.. seealso:: :ref:`flux units` 

..  _`flux units`:

What are the flux units?
------------------------

Flux units depend on the units of time and concentrations (of biomass and metabolites) provided in the input 
data file. For instance, if biomass units are in grams of cell dry weight by liter (gDW/L), metabolite concentrations are in millimolar (mM) and time is 
in hours (h), the estimated fluxes will be in mmol/gDW/h. Units should thus be carefully selected, and calculated fluxes must be interpreted consistently with the concentration units.

.. seealso:: :ref:`conc units` 

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

    * the plots of experimental vs simulated data for the best fit, which should be as close as possible.
    * the χ² statistical test results given in the log file (see below for help on interpreting the results).

.. seealso:: :ref:`chi2 test` and :ref:`bad fit` 

..  _`chi2 test`:

What is a χ² test?
------------------------------------------------------------------

A χ² test describes how well a model fits a set of observations. Measures of goodness of fit typically summarize the discrepancy between observed values and the values expected under the model used in PhysioFit (see the :doc:`models` section). It is calculated as the sum of differences between measured and simulated values, each squared and divided by the simulated value.
A good fit corresponds to small differences between measured and simulated values, thereby the χ² value is low. In contrast, a bad fit corresponds to large differences between simulations and measurements, and the χ² value is high. 

The resulting χ² value can then be compared with a χ² distribution to determine the goodness of fit. The p-value of one-tail χ² test is calculated by PhysioFit from the best fit and is given in the log file (have a look to the :doc:`usage` section). A p-value close to 0 means poor fitting, and a p-value close to 1 means good fitting (keeping in mind that a p-value very close to 1 can be an evidence that standard deviations might be overestimated). A 
p-value between 0.95 and 1 means the model fits the data good enough with respect to the standard deviations provided (at a 95% confidence level). PhysioFit provides an explicit meassage stating wether the flux data are satisfactorily fitted or not (at a 95% confidence interval).

..  _`bad fit`:

My data hasn't been correctly fitted. Why?
------------------------------------------------------------------

A possible reason to explain a bad fit is that standard deviations on measurements (concentration biomass and metabolites) is under-estimated, thereby making the χ² test too stringent. In this case, plots of measured and fitted data should be in agreement. Reliable estimated of standard deviation on measurements must be provided to PhysioFit (have a look to the :doc:`usage` section to see how to check and adjust this parameter).

Another possible reason to explain a bad fit is that a key asumption of the flux calculation method is not respected. Typically, cells might not be strictly in metabolic steady-state, i.e. with constant fluxes during the whole experiment. If this key asumption does not occur (e.g. cells are continuously adapting to their environment and fluxes change over time), PhysioFit will not be able to fit the data satisfactorily. In this case, evaluate wether the deviation is significant or not (e.g. based on the detailed χ² statistics or on the plot of fitted vs measured data), and evaluate the potential biases that would be introduced by interpreting (or not) these flux values.

In rare situations, it may also be because some parameters have to be tweaked to
help PhysioFit fitting the measurements, which results in obviously aberrant fits (e.g. with flat time-course profiles for all metabolites). This might happen for instance if some measurements are provided in units far from unity (e.g. 1e-5 M instead of 10 µM). If this situation happens, we suggest modifying the initial values of fluxes, or changing the units of input data, and re-run the flux calculation. For more info on the run parameters and how they may affect the fitting process,
please refer to section :ref:`physiofit parameters`.

If you believe the problem is in PhysioFit, we would greatly appreciate 
if you could open a new issue on our `issue tracker  <https://github.com/MetaSys-LISBP/PhysioFit/issues>`_.
   
I cannot start PhysioFit graphical user interface, can you help me?
-------------------------------------------------------------------

If you  installed PhysioFit following our standard procedure and that you are unable
to start PhysioFit by opening a terminal and typing :samp:`physiofit`, then there is indeed
something wrong. Do not panic, we are here to help!
Please follow this simple procedure:

1. The first step of the debugging process will be to get a *traceback*, i.e.
   a message telling us what is actually going wrong. You should see this message in the terminal you opened.

2. Read the traceback and try to understand what is going wrong:

   * If it is related to your system or your Python installation, you will need to ask some
     help from your local system administrator or your IT department so they could
     guide you toward a clean installation. Tell them that you wanted "to use the graphical
     user interface of PhysioFit, a Python 3.6 software" and what you did so far (installation),
     give them the traceback and a link toward the documentation. They should know what to do.
   * If you believe the problem is in PhysioFit or that your local system administrator
     told you so, then you probably have found a bug! We would greatly appreciate
     if you could open a new issue on our `issue tracker  <https://github.com/MetaSys-LISBP/PhysioFit/issues>`_.
     
I would like a new feature.
------------------------------------------------------------------

We would be glad to improve PhysioFit. Please get in touch with us `here 
<https://github.com/MetaSys-LISBP/PhysioFit/issues>`_ so we could discuss your problem.
