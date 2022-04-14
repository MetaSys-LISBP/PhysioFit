"""
Module containing the Command-Line Interface control logic.
"""

"""
###########
PROCEDURE #
###########

-> Initialize the IoHandler in local mode or galaxy mode.
-> Load data in. Two cases:
    - Data is tsv format
    - Data is config file (json) --> only in local mode?
In the json case, data must be imported using path_to_data parameter. Start by testing it's validity.
The loading of the data must be done using the galaxy_in function (to be written) if the galaxy argument is used. 
CLI will mostly only be utilized through Galaxy. 

-> Initialize the fitter
-> Run optimization
-> Run tests
-> Run monte-carlo if asked

######
CODE #
######

1. Write function to parse_args.

For this, we need to write out an exhaustive list of all the arguments. 
    
Optionals:
    DEVELOPPER
        - Galaxy logic (--galaxy, -g): Is the CLI being used on the galaxy platform.
        
    SECTION Data
        - Tsv file
        - Json file
    
    SECTION Basic parameters:
        - Lag phase (--lag, -l): Should lag phase be estimated. Add for True
        - Degradation (--deg, -d): Should degradation constants be taken into consideration. Add for True. Constants 
        should then be given in dict format (need specific function to parse this in command-line)
        - Monte Carlo (--montecarlo, -mc): Should sensitivity analysis be performed. Add for True. Number of iterations 
        should be given as an int.
?? How should path be handled??
    
    SECTION Advanced parameters:
        - Initial value (--vini, -i): Select an initial value for fluxes to estimate. Add and give value 
        (floating-point).
        - Weights (--weight, -w): Standard deviation on measurements. Add and give weights in dictionary format (use 
        same parsing function as for degradation).
        - Initial metabolite conc bounds (--conc_met_bounds, -cm): Bounds on initial metabolite concentrations. Add and 
        give bounds as a tuple.
        - Initial met flux bounds (--flux_met_bounds, -fm): Bounds on initial metabolite fluxes. Add and give bounds as
        a tuple.
        - Initial biom conc bounds (--conc_biom_bounds, -cb): Bounds on initial biomass concentrations. Add and give 
        bounds as a tuple.
        - Initial biom flux bounds (--flux_biom_bounds, -fb): Bounds on initial biomass fluxes. Add and give bounds as 
        a tuple.
        - Debug mode (--verbose, -v): Activate the debug logs. Add fo True.
        
2. Cli class

    A. Initialize arguments: home, run_home
Here we initialize the CLI object. Arguments are parsed and given to an argument parser attribute.



"""