"""
Module containing the Command-Line Interface control logic.
"""
import argparse
from ast import literal_eval
from copy import copy

from physiofit.base.io import IoHandler


def parse_args():
    """
    Parse arguments from user input.
    :return: Argument Parser object
    :rtype: class: argparse.ArgumentParser
    """

    parser = argparse.ArgumentParser(
        "Physiofit: Extracellular flux estimation software")

    # Parse data arguments (tsv + json)
    parser.add_argument("-t", "--data", type=str,
                        help="Path to data file in tabulated format (txt or tsv)")
    parser.add_argument("-c", "--config", type=str,
                        help="Path to config file in json format")

    # Parse basic parameters
    parser.add_argument("-l", "--lag", action="store_true",
                        help="Should lag phase be estimated. Add for True")
    parser.add_argument("-d", "--deg", type=str,
                        help="Should degradation constants be taken into "
                             "consideration. Add for True. Constants should"
                             " then be given in dictionary format")
    parser.add_argument("-mc", "--montecarlo", type=int,
                        help="Should sensitivity analysis be performed. Number of "
                             "iterations should be given as an integer"
                        )

    # Parse advanced parameters
    parser.add_argument("-i", "--vini", type=float,
                        help="Select an initial value for fluxes to estimate")
    parser.add_argument("-s", "--sd", type=str,
                        help="Standard deviation on measurements. Give weights in "
                             "dictionary format")
    parser.add_argument("-cm", "--conc_met_bounds", type=str,
                        help="Bounds on initial metabolite concentrations. "
                             "Give bounds as a tuple."
                        )
    parser.add_argument("-fm", "--flux_met_bounds", type=str,
                        help="Bounds on initial metabolite fluxes. "
                             "Give bounds as a tuple."
                        )
    parser.add_argument("-cb", "--conc_biom_bounds", type=str,
                        help="Bounds on initial biomass concentrations."
                             "Give bounds as a tuple."
                        )
    parser.add_argument("-fb", "--flux_biom_bounds", type=str,
                        help="Bounds on initial biomass fluxes. "
                             "Give bounds as a tuple."
                        )

    # Parse developer arguments
    parser.add_argument("-g", "--galaxy", action="store_true",
                        help="Is the CLI being used on the galaxy platform"
                        )
    parser.add_argument("-v", "--debug_mode",action="store_true",
                        help="Activate the debug logs"
                        )

    # Parse selective output path arguments (for galaxy implementation mostly)
    parser.add_argument(
        "-op", "--output_pdf", type=str,
        help="Path to output the pdf file containing plots"
    )
    parser.add_argument(
        "-of", "--output_fluxes", type=str,
        help="Path to output the flux results"
    )
    parser.add_argument(
        "-os", "--output_stats", type=str,
        help="Path to output the khi² test"
    )
    parser.add_argument(
        "-oc", "--output_config", type=str,
        help="Path to output the json config file"
    )


    return parser


class Cli:
    """
    Class to handle the io_handler initialization and CLI logic
    """

    def __init__(self):

        self.parser = parse_args()
        self.args = self.parser.parse_args()
        self.io_handler = None

        if not self.args.config or self.args.galaxy:
            self.fitter_args = self._initialize_fitter_args()

    def _initialize_fitter_args(self) -> dict:
        """
        Setup arguments that shall be fed to the fitter object

        :return: fitter arguments
        """

        fitter_args = copy(vars(self.args))

        # Remove args that where not passed in by the cli
        for key in self.args.__dict__.keys():
            if self.args.__dict__[key] is None:
                del fitter_args[key]

        # Clean up the arguments and give them fitter parameter names
        if "lag" in fitter_args.keys():
            fitter_args["t_lag"] = fitter_args.pop("lag")
        if "montecarlo" in fitter_args.keys() \
                and "config" not in fitter_args.keys():
            fitter_args["mc"] = True
            fitter_args["iterations"] = fitter_args.pop("montecarlo")

        # Transform the strings into the right types
        for arg in ["deg", "conc_met_bounds", "flux_met_bounds",
                    "conc_biom_bounds", "flux_biom_bounds"]:
            if arg in fitter_args.keys():
                fitter_args[arg] = literal_eval(fitter_args.pop(arg))

        # Remove arguments that are data or config files or output paths
        for arg in [
            "data", "config", "galaxy", "output_pdf",
            "output_plots", "output_fluxes", "output_stats",
            "output_config"
        ]:
            if arg in fitter_args.keys():
                del fitter_args[arg]

        return fitter_args

    def start_cli(self):
        """
        Launch the Command Line Interface in local or Galaxy mode
        :return: None
        """

        print("Starting CLI")
        if not self.args.galaxy:
            print("CLI started in local mode")
            self.local_launch()
        else:
            print("CLI started in Galaxy mode")
            self.galaxy_launch()

    def local_launch(self):
        """
        Logic to handle the Command Line Interface in local mode
        :return: None
        """

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
        """
        Logic to handle the Command Line Interface in Galaxy mode
        :return: None
        """

        self.io_handler = IoHandler("galaxy")

        # Ensure data is given through the galaxy UI with the config file
        if self.args.config and not self.args.data:
            raise DataError(
                "To run on Galaxy using a config file, "
                "Physiofit 2.0 needs the data to be given as well."
            )
        elif self.args.config and self.args.data:
            print("Json config file detected")
            self.io_handler.galaxy_json_launch(self.args.config,
                                               self.args.data)
        # Handle case where only data is given
        elif self.args.data and not self.args.config:
            print("Data file detected")
            self.io_handler.galaxy_in(self.args.data, **self.fitter_args)
        else:
            raise DataError("No config file or data file detected")

        print("Fitter initialized")
        self.run()

    def run(self):
        """
        Launch optimization, monte carlo analysis, stat test and output results
        :return: None
        """

        print("Running optimization...")
        self.io_handler.fitter.optimize()

        if self.io_handler.fitter.mc:
            print("Running sensitivity analysis...")
            self.io_handler.fitter.monte_carlo_analysis()

        print("Running khi2 test...")
        self.io_handler.fitter.khi2_test()

        print("Exporting results...")
        if self.io_handler.input_source == "galaxy":
            self.io_handler.output_pdf(self.args.output_pdf)
            self.io_handler.output_report(
                [self.args.output_stats, self.args.output_fluxes]
            )
            self.io_handler._generate_run_config(self.args.output_config)

        else:
            self.io_handler.local_out("data", "plot", "pdf")
        print("Done!")


class DataError(Exception):
    pass


def main():
    """
    Main function to call CLI and launch CLI
    :return: None
    """
    physiofit_cli = Cli()
    physiofit_cli.start_cli()
