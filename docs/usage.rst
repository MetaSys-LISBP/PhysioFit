Usage
=====

Required input data file
------------------------

The input data must be a tabulated file (.tsv extension) to be imported into the application. The structure of
the file is as follows:

==== ===== ======
Time   X    Met1
==== ===== ======
 1    0.2     x
 2    0.3     x
 3    0.4     x
etc   ...    ...
==== ===== ======

The time and biomass(X) columns are mandatory, as is at least one metabolite column. If the biomass and metabolite
concentrations were sampled at different moments, you can still group them together in the same table.

The json config file
---------------------

The json config file presents an organisation similar to a python dictionary. You can find an example of one `here <>`

PhysioFit Parameters
--------------------

Here is a list of the different parameters that can be modified when setting up your run:

    * **Basic parameters:**
        - *Lag*: Should the length of the lag phase be estimated during optimization
        - *Degradation*: Should degradation constants be used during optimization (if yes then they should be given in
          dictionary format).
        - *Sensitivity analysis (Monte Carlo)*: Should the sensitivity analysis be performed on optimized parameters (if
          true then input the number of iterations).

    * **Advanced parameters:**
        - *vini*: Initial value for parameters to estimate. This will be the starting point for the optimization.
        - *Weights*: Standard deviation to apply during cost calculation (see eq. 10 in :doc:`method`)
        - *Metabolite concentration bounds*: Bounds to apply to the metabolite concentrations during optimization to
          reduce the range of possibilities.
        - *Flux bounds*: Bounds to apply on fluxes during optimization to reduce the range of possibilities
        - *Biomass initial concentration bounds*: Bounds to apply on the starting concentration for biomass to reduce
          the range of possibilities.
        - *Growth rate bounds*: Bounds to apply on the growth rate to reduce the range of possibilities.

Result files
-------------