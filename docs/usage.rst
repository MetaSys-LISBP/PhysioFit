Tutorial
========

Required input data file
------------------------

The input data must be a tabulated file (.tsv extension) with the following structure:

+------+-----+---------+-----+
| time |  X  | Glucose | ... |
+======+=====+=========+=====+
| 1    | 0.2 |   13    | ... |
+------+-----+---------+-----+
| 2    | 0.3 |   12    | ... |
+------+-----+---------+-----+
| 2.5  | 0.4 |   10    | ... |
+------+-----+---------+-----+
| ...  | ... |   ...   | ... |
+------+-----+---------+-----+

Columns **time** and **X** (biomass concentration) are mandatory, as is at least one metabolite column (**Glucose** in this example). If the biomass and metabolite
concentrations were sampled at different moments, you can still group them together in the same table (with an empty value when no data is available).

.. note:: Flux units depend on the units of time and concentrations (of biomass and metabolites) provided in the input
             data file. For instance, if biomass units are in grams of cell dry weight by liter (gDW/L), metabolite concentrations are in millimolar (mM) and time is
             in hours (h), the estimated fluxes will be in mmol/gDW/h. Units should thus be carefully selected, and calculated fluxes must be interpreted consistently with the concentration units.

.. warning:: To limit any numerical instabilities, provide values in a range not too far from unity (e.g. if a metabolite
             concentration is 2 mM, provide the value directly in mM and not as 2e-3 M). The concentration of different metabolites can
             be provided using different units, but a single unit must be used for a given metabolite.

.. _json_config
The json configuration file
---------------------

The json configuration file contains all parameters required to calculate fluxes for a given experiment, and is 
thus key to ensure **reproducibility** of the flux calculation process. A configuration file is generated automatically by 
PhysioFit during flux calculation; it can also be created or edited manually.

Its structure is similar to a python dictionary. You can find an example of configuration file `here
<https://github.com/MetaSys-LISBP/PhysioFit/blob/dev_v2.0/config_example_file.json>`_.

For a description of all calculation parameters, check the section below.

.. _PhysioFit parameters:

Flux calculation parameters
--------------------

This section presents the different parameters that can be adapted to match your own experiments. Additional information on each parameter can be found in the :doc:`method` section.

    * **Basic parameters**:
        - **Lag**: consider (or not) a lag phase (i.e. without growth) during flux calculation (see eq. 1 in :doc:`method` section).
        - **Degradation**: (first-order) degradation constants to be used during flux calculation (constants should be given as a
          Python dictionary, see eq. 2 in :doc:`method` section).
        - **Sensitivity analysis (Monte Carlo)**: Estimate (or not) the precision on calculated fluxes, using a Monte Carlo sensitivity analysis. If
          selected, then provide the number of monte carlo iterations. A higher number of iterations will give more accurate confidence
          intervals on the estimated parameters, but will slow down calculations. The default number of
          iterations (100) is sufficient in most situations.

    * **Advanced parameters**:
        - **Initial flux values (qM and µ)**: Initial value for fluxes and growth rate to estimate (see eqs. 1 and 2 in :doc:`method` section). Default: 0.2
        - **Standard deviation**: Standard deviation on measurements (see weight in eq. 10 of the :doc:`method` section). A smaller value
          will augment the cost of the corresponding data during the optimization, thereby forcing an improvement of the fit accuracy for this data, but degrading the fit accuracy for the other measurements. Defaults: 0.02 for biomass, and 0.5 for metabolites.
        - **Bounds on initial metabolite concentrations (Mi0)**: Minimal and maximal values of the initial concentration of metabolites. Default: [1e-06, 50]
        - **Flux bounds**: Minimal and maximal values of fluxes. Default:
          [-50, 50]
        - **Bounds on initial biomass concentration (X0)**: Minimal and maximal values of the initial concentration of biomass. Default: [1e-03, 2]
        - **Bounds on growth rate (µ)**: Minimal and maximal values of the growth rate. Default: [1e-4, 50]
        - **Verbose logs**: Should debug information be written in log file. Useful in case of trouble (please join it to the issue on github). Default: False

.. note:: Initial values and bounds should be carefully chosen. Ideally, initial values should be in the range of values used in the experiment. Well-defined bounds will enhance robustness and speed of the flux calculation process. The default
          bounds are sufficient in most cases, but may still be defined by the user when needed (e.g. the higher bound on initial metabolite concentrations should be increased if the initial concentration of substrate is higher than 50, since it is the maximal value allowed by default.).

Output files
---------------------------

The following files are generated by PhysioFit in the output directory:

    * **config_file.json**: This is the configuration file containing all parameters used for the last run.
    * **flux_results.tsv**: Flux calculation results, i.e. fluxes and initial metabolite concentrations for the best fit, with associated precision.
    * **stat_results.tsv**: Results from the χ² statistical test.
    * **log.txt**: The run log containing information on how the run went.
    * **plots.pdf**: pdf file containing plots of simulated and measured data.
    * A number of **.svg files** containing the individual plots of simulated and measured data.

.. note::
   * When the data file (and not a json configuration file) is directly used as input in the Graphical User Interface, PhysioFit cannot get the path directly from the file metadata. For this reason, *path_to_data* is set to None in the generated *config_file.json*.

   * PhysioFit silently overwrites (results and log) files if they already exist. So take care to copy your results elsewhere if you want to protect them from overwriting.

The quality of the fit must be checked before interpreting the estimated fluxes by checking the *test_results.tsv* 
file, which contains the detailed χ² statistical test results and a clear status on the quality of the fit (based on a 95% confidence interval). The generated plots also help to visualize how accurately the simulated data fits the
experimental measurements. Finally, the confidence intervals estimated using the Monte-Carlo approach provides quantitative information on the precision of the estimated fluxes (mean, median, standard deviation and 95% confidence interval). 
Have a look to the :doc:`faq` section for help on interpreting the statistical results.
