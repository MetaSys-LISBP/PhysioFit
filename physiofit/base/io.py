"""Module to handle inputs and outputs for PhysioFit"""
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from pandas import DataFrame, read_csv

from physiofit.base.fitter import PhysioFitter

mod_logger = logging.getLogger("PhysioFit.base.io")


class IoHandler:
    """
    Input/Output class that handles the former and initializes the PhysioFitter component object. It is the
    preferred interface for interacting with the PhysioFit package.

    :param source: Input source
    """

    allowed_keys = {
        "vini", "conc_biom_bounds", "flux_biom_bounds", "conc_met_bounds", "flux_met_bounds", "weight", "t_lag", "deg",
        "iterations", "mc", "debug_mode"
    }

    def __init__(self, source='local'):

        self.input_source = source
        self._output_type = []
        self._figures = []
        self.fitter = None
        self.output = None
        self.data = None
        self.names = None
        self.simulated_data_sds = None
        self.simulated_data = None
        self.experimental_data = None
        self.home_path = None
        self.data_path = None
        self.res_path = None
        self.has_config_been_read = False

    @staticmethod
    def _read_data(path_to_data: str) -> DataFrame:
        """
        Read initial data file (csv or tsv)

        :param path_to_data: str containing the relative or absolute path to the data
        :return: pandas DataFrame containing the data
        """

        data_path = Path(path_to_data).resolve()

        if data_path.suffix == ".tsv":
            data = read_csv(str(data_path), sep="\t")
        elif data_path.suffix == ".csv":
            data = read_csv(str(data_path), sep=";")
        else:
            if not data_path.exists():
                raise ValueError(f"{data_path} is not a valid file")
            else:
                raise TypeError(f"{data_path} is not a valid format. Accepted formats are .csv or .tsv")

        IoHandler._verify_data(data)
        return data

    @staticmethod
    def _verify_data(data: DataFrame):
        """
        Perform checks on DataFrame returned by the _read_data function

        :param data: pandas DataFrame containing the data
        :return: None
        """

        if not isinstance(data, DataFrame):
            raise TypeError("There was an error reading the data: DataFrame has not been generated")

        for col in ["time", "X"]:
            if col not in data.columns:
                raise ValueError(f"Column {col} is missing from the dataset")

        if len(data.columns) <= 2:
            raise ValueError(f"Data does not contain any metabolite columns")

        for col in data.columns:
            if data[col].dtypes != np.int64 and data[col].dtypes != np.float64:
                raise ValueError(f"Column {col} has values that are not of numeric type")

    @staticmethod
    def generate_config_file(destination_path: str):
        """
        Generate the default configuration file

        :param destination_path: path to configuration file
        :return: None
        """

        config = {
            key: ""
            for key in IoHandler.allowed_keys
        }

        config["vini"] = None
        config["mc"] = None
        config["iterations"] = None
        config["conc_biom_bounds"] = None
        config["flux_biom_bounds"] = None
        config["conc_met_bounds"] = None
        config["flux_met_bounds"] = None
        config["weight"] = None
        config["deg"] = None
        config["t_lag"] = None
        config["debug_mode"] = False
        config.update(
            {"path_to_data": None}
        )

        dest_path = Path(destination_path) / "config_file.json"

        with open(str(dest_path), "w") as conf:
            json.dump(config, conf, indent=4, sort_keys=False)

    def local_in(self, data_path: str, **kwargs: dict):
        """
        Function for reading data and initializing the fitter object in local context

        :param data_path: path to data
        :param kwargs: supplementary keyword arguments are passed on to the PhysioFitter object
        :return: None
        """

        # Store the data path
        data_path = Path(data_path).resolve()
        if not data_path.is_absolute():
            data_path = data_path.absolute()

        # Initialize data and set up destination directory
        if self.data is not None:
            raise ValueError(f"It seems data is already loaded in. Data: \n{self.data}\nHome path: {self.home_path}")
        elif self.input_source == "local":
            self.home_path = data_path.parent.resolve()
            self.res_path = self.home_path / (self.home_path.name + "_res")
            if not self.res_path.is_dir():
                self.res_path.mkdir()
            if not self.home_path.exists():
                raise KeyError(f"Input file does not exist. Path: {self.home_path}")
            self.data = IoHandler._read_data(data_path)
            self.data = self.data.sort_values("time", ignore_index=True)
            self.names = self.data.columns[1:].to_list()
        else:
            raise IOError(f"Wrong input source selected. Source: {self.input_source}")

        self.initialize_fitter(kwargs)

        if not self.has_config_been_read:
            self._generate_run_config()

    def _generate_run_config(self):
        """
        Generate configuration file from parameters of the last run

        :return: None
        """

        to_dump = {}
        for key, value in self.fitter.__dict__.items():
            if key in self.allowed_keys:
                if isinstance(value, np.ndarray):
                    to_dump.update({key: value.tolist()})
                else:
                    to_dump.update({key: value})

        to_dump.update(
            {"path_to_data": str(self.data_path)}
        )

        with open(str(self.res_path / "config_file.json"), "w") as conf:
            json.dump(to_dump, conf, indent=4, sort_keys=True)
        self.fitter.logger.info(f"\nConfiguration file saved in: {self.res_path / 'config_file.json'}")

    def launch_from_json(self, json_file: str or bytes):
        """
        Launch the run using a json file as input

        :param json_file: json file containing run parameters. Can be json string or file-like or path to file
        :return: None
        """

        config = self.read_json_config(json_file)
        self.has_config_been_read = True

        # Get the data path and remove from config dict to ensure no wrong key errors are returned during fitter
        # initialization
        data_path = config["path_to_data"]
        del config["path_to_data"]

        self.local_in(data_path, **config)

    @staticmethod
    def read_json_config(json_file: str or bytes) -> dict:
        """
        Import json configuration file and parse keyword arguments

        :param json_file: path to the json file or json file
        :return config: Dictionnary containing arguments parsed from json file
        """

        # Load config file
        if type(json_file) == str:
            try:
                path = Path(json_file).resolve()
                json_file = open(str(path))
                config = json.load(json_file)
            except OSError:
                config = json.loads(json_file)
        else:
            config = json.load(json_file)

        # Convert lists to tuples for the bounds
        if type(config["conc_biom_bounds"]) == list:
            config["conc_biom_bounds"] = tuple(config["conc_biom_bounds"])
        if type(config["conc_met_bounds"]) == list:
            config["conc_met_bounds"] = tuple(config["conc_met_bounds"])
        if type(config["flux_met_bounds"]) == list:
            config["flux_met_bounds"] = tuple(config["flux_met_bounds"])
        if type(config["flux_biom_bounds"]) == list:
            config["flux_biom_bounds"] = tuple(config["flux_biom_bounds"])

        # Remove None values from config dict so that defaults are used on fitter init
        keys_to_del = [key for key, value in config.items() if value is None]
        for key in keys_to_del:
            del config[key]

        return config

    def local_out(self, *args):
        """
        Function for coordinating exports in local context

        :param args: type of output to generate (data, plot and/or pdf)
        """

        for arg in args:
            if arg not in ["data", "plot", "pdf"]:
                raise ValueError(f"Detected argument is not an output type: {arg}.\n"
                                 f"Accepted output types are: data, plot, pdf")
            else:
                self._output_type.append(arg)

        if "data" in self._output_type:
            self._output_report()

        if "plot" in self._output_type:
            self.plot_data()
            self._output_plots()

        if "pdf" in self._output_type:
            self._output_pdf()

        self._generate_run_config()

    def initialize_fitter(self, kwargs: dict = None):
        """
        Initialize a PhysioFitter object

        :param kwargs: Keyword arguments for fitter initialization
        :return: None
        """

        wrong_keys = []

        # Initialize fitter
        self.fitter = PhysioFitter(self.data)

        # Initialize fitter logger
        file_handle = logging.FileHandler(self.res_path / "log.txt", "w+")
        try:
            if kwargs["debug_mode"]:
                file_handle.setLevel(logging.DEBUG)
            else:
                file_handle.setLevel(logging.INFO)
        except KeyError:
            file_handle.setLevel(logging.INFO)
        except Exception:
            raise "An error has occurred while initializing the log level"

        file_handle.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        if not self.fitter.logger.handlers:
            self.fitter.logger.addHandler(file_handle)
        self.fitter.logger.debug(f"Logger handlers: {self.fitter.logger.handlers}")

        if kwargs:
            # Initialize fitter parameters from kwargs
            for key, value in kwargs.items():
                self.fitter.logger.debug(f"Key: {key}, value: {value}\n")
                if key in IoHandler.allowed_keys:
                    self.fitter.__dict__.update({key: value})
                else:
                    wrong_keys.append(key)
        self._initialize_fitter_vars()
        self.fitter.verify_attrs()
        self.fitter.logger.debug(f"Fitter attribute dictionary:\n{self.fitter.__dict__}")
        if wrong_keys:
            raise KeyError(f"Some keyword arguments were not valid: {wrong_keys}")

    def _initialize_fitter_vars(self):
        """
        Initialize fitter variables
        :return: None
        """

        self.fitter.initialize_vectors()
        self.fitter.initialize_weight_matrix()
        self.fitter.initialize_bounds()
        self.fitter.initialize_equation()

    def _output_pdf(self):
        """Handle the creation and output of a pdf file containing fit results as a plot"""

        if not self.home_path:
            raise RuntimeError("No home path detected. Was data loaded in correctly?")
        if not self._figures:
            raise RuntimeError("Plots have not been created. Please launch the plot_data() function first")

        try:
            with PdfPages(rf"{self.res_path}\plots.pdf") as pdf:
                for fig in self._figures:
                    pdf.savefig(fig[1])
        except Exception:
            raise "Error while generating the pdf file"

    def _output_plots(self):
        """Handle the creation and export of the different plots in svg format"""

        if not self.home_path:
            raise RuntimeError("No home path detected. Was data loaded in correctly?")
        if not self._figures:
            raise RuntimeError("Plots have not been created. Please launch the plot_data() function first")

        try:
            for fig in self._figures:
                fig[1].savefig(rf"{self.res_path}\{fig[0]}.svg")
        except Exception:
            raise RuntimeError("Unknown error while generating output")

    def _output_report(self):
        """
        Handle creation and export of the report containing stats from monte carlo analysis of optimization
        parameters
        """

        self.fitter.logger.debug(f"Parameter Stats:\n{self.fitter.parameter_stats}")

        # Get optimization results as dataframe
        opt_data = DataFrame.from_dict(self.fitter.parameter_stats,
                                       orient="columns")

        # Use IDs to clarify which parameter is described on each line
        opt_data.index = self.fitter.ids
        opt_data.to_csv(fr"{self.res_path}\flux_results.tsv", sep="\t")

        if isinstance(self.fitter.khi2_res, DataFrame):
            with open(fr"{self.res_path}\stat_results.tsv", "w+") as stat_out:
                stat_out.write("==================\n"
                               "Khi² test results\n"
                               "==================\n\n")
                stat_out.write(self.fitter.khi2_res.to_string(header=False, justify="center"))
                if self.fitter.khi2_res.at["p_val", "Values"] < 0.95:
                    stat_out.write(
                        f"\n\nAt level of 95% confidence, the model fits the data good enough with respect to the "
                        f"provided measurement SD. Value: "
                        f"{self.fitter.khi2_res.at['p_val', 'Values']}")

                else:
                    stat_out.write(
                        f"\n\nAt level of 95% confidence, the model does not fit the data good enough with respect to "
                        f"the provided measurement SD. "
                        f"Value: {self.fitter.khi2_res.at['p_val', 'Values']}\n")

    def _get_plot_data(self):
        """
        Prepare data for plotting
        """

        # Initialize vectors and data for plotting
        if self.fitter.time_vector is not None:
            x = self.fitter.time_vector
        else:
            raise ValueError("PhysioFitter time_vector has not been initialized. "
                             "Have you loaded in the data correctly?")
        if self.fitter.experimental_matrix is not None:
            exp_mat = self.fitter.experimental_matrix
        else:
            raise ValueError("PhysioFitter object does not have experimental data loaded in")
        if self.fitter.simulated_matrix is not None:
            sim_mat = self.fitter.simulated_matrix
        else:
            raise ValueError("PhysioFitter simulated data does not exist yet")
        if self.fitter.matrices_ci is not None:
            sim_mat_ci = self.fitter.matrices_ci
            self.simulated_data_low_ci = DataFrame(columns=self.names, index=x, data=sim_mat_ci["lower_ci"])
            self.simulated_data_high_ci = DataFrame(columns=self.names, index=x, data=sim_mat_ci["higher_ci"])
        else:
            self.fitter.logger.warning(
                "Monte Carlo analysis has not been done, confidence intervals will not be computed")

        self.experimental_data = DataFrame(columns=self.names, index=x, data=exp_mat)
        self.simulated_data = DataFrame(columns=self.names, index=x, data=sim_mat)

    def plot_data(self, display=False):
        """
        Plot the data

        :param display: should plots be displayed
        """

        self._get_plot_data()
        self._draw_plots(display)

    def _draw_plots(self, display: bool):
        """
        Draw plots and assign them to the _figures attribute for later access in pdf and plot generation functions

        :param display: Should plots be displayed or not on creation
        """

        for element in self.names:

            # Initialize figure object
            fig, ax = plt.subplots()
            fig.set_size_inches(9, 5)

            # Draw experimental points onto the Axe
            exp = ax.scatter(self.experimental_data.index,
                             self.experimental_data[element],
                             marker='o', color="orange")
            exp.set_label(f"Exp {element}")

            # Draw the simulated line onto the Axe
            sim_line, = ax.plot(self.simulated_data.index,
                                self.simulated_data[element],
                                linestyle='--')
            sim_line.set_label(f"Sim {element}")

            # If monte carlo analysis has been done add the corresponding sd to the plot as a red area
            if self.fitter.matrices_ci is not None:
                ax = self._add_sd_area(element, ax)

            # Finishing touches
            ax.set(xlim=0, ylim=0, xlabel="Time", ylabel="Concentration")
            ax.legend()
            ax.set_title(f"{element}")
            fig.tight_layout()

            # Add the figure with the metabolite name as index to the _figures attribute for later use
            self._figures.append((element, fig))

        if display:
            for _ in self._figures:
                plt.show()

    def _add_sd_area(self, element: str, ax: plt.Axes):
        """
        Draw red area around the fitting line to show confidence intervals

        :param element: Which variable is being plotted
        :param ax: axes on which to draw the area before returning to mother function
        """

        y1 = self.simulated_data_low_ci[element].to_numpy()
        y2 = self.simulated_data_high_ci[element].to_numpy()
        x = self.simulated_data.index
        ax.fill_between(x, y1, y2, alpha=.3, linewidth=0, color="red")
        return ax
