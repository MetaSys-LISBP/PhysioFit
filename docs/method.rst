Method
======

Initial concentrations of species (i.e. concentrations of biomass - X - and metabolites - M_i - at t=0) and fluxes
(exchange fluxes - q_(M_i ) - and growth rate - µ -) are estimated by fitting to time-course measurements of
experimental data. Simulations, i.e. calculation of X(t) and M_i (t), are performed using analytical functions, as
detailed below.

The general model, taking into accounts both degradation and lag phase, can be described using the following system of
ordinary differential equations:

.. image:: _static/equations/eq1.png

.. image:: _static/equations/eq2.png

with qMi being positive (negative) when M_i is produced (consumed). The sign of q_(M_i ) can thus be used to
automatically identify products and substrates in high throughput workflows for automated functional analysis of
metabolic systems.

Integrating equations 1-2 provides the following analytical functions:

.. image:: _static/equations/eq3.png

.. image:: _static/equations/eq4.png

In the absence of a lag phase (i.e. t_lag=0), equations 3-4 simplifies to:

.. image:: _static/equations/eq5.png

.. image:: _static/equations/eq6.png

In the absence of degradation (i.e. k = 0), eq. 4 simplifies to:

.. image:: _static/equations/eq7.png

In the absence of both degradation and lag (i.e. t_lag=0 and k=0), equations 3-4 simplifies to:

.. image:: _static/equations/eq8.png

.. image:: _static/equations/eq9.png

Parameter estimation is performed using Scipy's L-BFGS-B method (see
`scipy.optimize.minimize <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`_ for more
information) by minimizing the following cost function:

.. image:: _static/equations/eq10.png

where *sim* is the simulated data, *meas* denotes measurements, and *weight* is a weighting factor (e.g., standard
deviation on measurements).

A global sensitivity analysis (Monte-Carlo approach) is implemented to evaluate the precision of the estimated
parameters, plots are generated for visual inspection of the fitting quality, and a khi2 test is performed to assess the
statistical goodness of fit.

Finally, PhysioFit includes routines to estimate the first-order degradation constants from time-course metabolite
concentrations measured in the medium without cells. In this situation, simulations are performed using the following
equation:

.. image:: _static/equations/eq11.png
