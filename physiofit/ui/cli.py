"""
Module containing the Command-Line Interface control logic.
"""
import argparse
from ast import literal_eval
from copy import copy

from physiofit.base.io import IoHandler

"""
###########
PROCEDURE #
###########

-> Initialize the IoHandler in local mode or galaxy mode.
-> Load data in. Two cases:
    - Data is tsv format
    - Data is config file (json)
Since there are two modes for the cli (galaxy and local), there are two ways of loading in the data:
    - In local mode, if a config file is given, the path_to_data parameter must lead to the tsv data file. 
    - In Galaxy mode, if a config file is given, the tsv file must also be given in.
In Galaxy mode, the parameters that are given in the config file must be entered into the rest of the GUI. 

-> Initialize the IoHandler. Two cases:
    - Source = local
Here we use a classical approach, with the IoHandler in local mode.  
    - Source = galaxy

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

"""


def parse_args():
    """
    Parse arguments from user input.
    :return: Argument Parser object
    :rtype: class: argparse.ArgumentParser
    """

    parser = argparse.ArgumentParser("Physiofit: Extracellular flux estimation software")

    # Parse data arguments (tsv + json)
    parser.add_argument("-t", "--data", type=str, help="Path to data file in tabulated format (txt or tsv)")
    parser.add_argument("-c", "--config", type=str, help="Path to config file in json format")

    # Parse basic parameters
    parser.add_argument("-l", "--lag", action="store_true", help="Should lag phase be estimated. Add for True")
    parser.add_argument("-d", "--deg", type=str, help="Should degradation constants be taken into "
                                                                 "consideration. Add for True. Constants should"
                                                                 " then be given in dictionary format")
    parser.add_argument("-mc", "--montecarlo", type=int, help="Should sensitivity analysis be performed. Number of"
                                                              "iterations should be given as an integer")

    # Parse advanced parameters
    parser.add_argument("-i", "--vini", type=float, help="Select an initial value for fluxes to estimate")
    parser.add_argument("-s", "--sd", type=str, help="Standard deviation on measurements. Give weights in "
                                                     "dictionary format")
    parser.add_argument("-cm", "--conc_met_bounds", type=str,
                        help="Bounds on initial metabolite concentrations. Give bounds as a tuple.")
    parser.add_argument("-fm", "--flux_met_bounds", type=str,
                        help="Bounds on initial metabolite fluxes. Give bounds as a tuple.")
    parser.add_argument("-cb", "--conc_biom_bounds", type=str,
                        help="Bounds on initial biomass concentrations. Give bounds as a tuple.")
    parser.add_argument("-fb", "--flux_biom_bounds", type=str,
                        help="Bounds on initial biomass fluxes. Give bounds as a tuple.")

    # Parse developer arguments
    parser.add_argument("-g", "--galaxy", action="store_true", help="Is the CLI being used on the galaxy platform")
    parser.add_argument("-v", "--verbose", help="Activate the debug logs")

    return parser


class Cli:

    def __init__(self):

        self.parser = parse_args()
        self.args = self.parser.parse_args()
        self.io_handler = None

        if not self.args.config or self.args.galaxy:
            self.fitter_args = self._initialize_fitter_args()

    def _initialize_fitter_args(self):

        fitter_args = copy(vars(self.args))

        # Remove args that where not passed in by the cli
        for key in self.args.__dict__.keys():
            if self.args.__dict__[key] is None:
                del fitter_args[key]

        # Clean up the arguments and give them fitter parameter names
        if "lag" in fitter_args.keys():
            fitter_args["t_lag"] = fitter_args.pop("lag")
        if "montecarlo" in fitter_args.keys() and "config" not in fitter_args.keys():
            fitter_args["mc"] = True
            fitter_args["iterations"] = fitter_args.pop("montecarlo")

        # Transform the strings into the right types
        for arg in ["deg", "conc_met_bounds", "flux_met_bounds", "conc_biom_bounds", "flux_biom_bounds"]:
            if arg in fitter_args.keys():
                fitter_args[arg] = literal_eval(fitter_args.pop(arg))

        # Remove arguments that are data or config files
        for arg in ["data", "config", "galaxy"]:
            if arg in fitter_args.keys():
                del fitter_args[arg]

        return fitter_args

    def start_cli(self):

        print("Starting CLI")
        if not self.args.galaxy:
            print("CLI started in local mode")
            self.local_launch()
        else:
            print("CLI started in Galaxy mode")
            self.galaxy_launch()

    def local_launch(self):

        self.io_handler = IoHandler("local")
        if self.args.config:
            print("Json config file detected")
            self.io_handler.launch_from_json(self.args.config)
        elif self.args.data:
            print("Data file detected")
            self.io_handler.local_in(self.args.data, **self.fitter_args)
        else:
            raise DataError("No config file or data file detected")
        print("Fitter initialized")
        self.run()

    def galaxy_launch(self):

        self.io_handler = IoHandler("galaxy")
        if self.args.config and not self.args.data:
            raise DataError("To run on Galaxy using a config file, Physiofit 2.0 needs the data to be given as well.")
        elif self.args.config and self.args.data:
            print("Json config file detected")
            self.io_handler.galaxy_json_launch(self.args.config, self.args.data)
        elif self.args.data and not self.args.config:
            print("Data file detected")
            self.io_handler.galaxy_in(self.args.data, **self.fitter_args)
        else:
            raise DataError("No config file or data file detected")
        print("Fitter initialized")
        self.run()

    def run(self):

        print("Running optimization...")
        self.io_handler.fitter.optimize()
        if self.io_handler.fitter.mc:
            print("Running sensitivity analysis...")
            self.io_handler.fitter.monte_carlo_analysis()
        print("Running khi2 test...")
        self.io_handler.fitter.khi2_test()
        print("Exporting results...")
        self.io_handler.local_out("data", "plot", "pdf")
        print("Done!")

class DataError(Exception):
    pass

def main():

    physiofit_cli = Cli()
    physiofit_cli.start_cli()
