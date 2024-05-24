"""Module to handle inputs and outputs for PhysioFit"""

from __future__ import annotations
import importlib
import logging
import os
import sys
from pathlib import Path
from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from pandas import DataFrame, read_csv, concat
import yaml

import physiofit.models.base_model
from physiofit import __file__
from physiofit.base.fitter import PhysioFitter
from physiofit.models.base_model import StandardDevs, Bounds, Model

# Switch matplotlib logger to higher level to not get debug logs in root logger
logging.getLogger("matplotlib").setLevel(logging.WARNING)

logger = logging.getLogger(f"physiofit.{__name__}")


class IoHandler:
    """
    Input/Output class that handles the former and initializes the
    PhysioFitter component object. It is the preferred interface for
    interacting with the PhysioFit package.

    """

    allowed_keys = {
        "sd", "model", "iterations", "mc", "debug_mode"
    }

    def __init__(self):

        self._output_type = []
        self.figures = []
        self.models = []
        self.fitter = None
        self.output = None
        self.data = None
        self.names = None
        self.simulated_data_sds = None
        self.simulated_data = None
        self.experimental_data = None
        self.wkdir = None
        self.data_path = None
        self.res_path = None
        self.has_config_been_read = False
        self.experiments = None
        self.multiple_experiments = False

    @staticmethod
    def read_data(data: str) -> DataFrame:
        """
        Read initial data file (csv or tsv)

        :param data: str containing the relative or
                             absolute path to the data
        :return: pandas DataFrame containing the data
        """
        try:
            if isinstance(data, str):
                data_path = Path(data).resolve()
                # .dat file type for galaxy implementation
                if data_path.suffix in [".txt", ".tsv", ".dat"]:
                    data = read_csv(str(data_path), sep="\t")
                elif data_path.suffix == ".csv":
                    data = read_csv(str(data_path), sep=";")
                else:
                    if not data_path.exists():
                        raise ValueError(f"{data_path} is not a valid file")
                    else:
                        raise TypeError(
                            f"{data_path} is not a valid format. "
                            f"Accepted formats are .csv, .txt or .tsv"
                        )
            elif issubclass(type(data), BytesIO):
                data = read_csv(data, sep="\t")
            else:
                raise TypeError(f"Input data file is not of right type. "
                                f"Accepted types: file-like (bytes) or string")
        except Exception:
            logger.exception("There was an error while reading the data")
            raise IOError(
                "Error while reading data. Please ensure you have the right "
                "file format (txt, tsv or bytes)"
            )

        IoHandler._verify_data(data)
        return data

    def select_model(
            self, name: str,
            data: pd.DataFrame = None
    ) -> Model:
        """
        Select a model from the list of models in the model folder of the
        package src directory
        """

        self.get_models(data)
        for model in self.models:
            if model.name == name:
                return model

    @staticmethod
    def read_model(
            model_file: str
    ):
        """
        Import and return the model class from .py file containing the model.

        .. warning: ONLY USE THIS FUNCTION ON TRUSTED PYTHON FILES. Reading
        code from untrusted sources can lead to propagation of viruses and
        compromised security.

        :param model_file: path to the model.py file to import
        """

        if not Path(model_file).is_file() or Path(model_file).suffix != ".py":
            raise ValueError(
                "The given path is not valid. The path must point to a .py "
                "file containing the module from which the model will be "
                "loaded."
            )

        spec = importlib.util.spec_from_file_location("module_to_import",
                                                      fr"{model_file}")
        module = importlib.util.module_from_spec(spec)
        sys.modules["module_to_import"] = module
        spec.loader.exec_module(module)
        model_class = getattr(module, "ChildModel")

        return model_class

    @staticmethod
    def _verify_data(
            data: DataFrame
    ):
        """
        Perform checks on DataFrame returned by the _read_data function

        :param data: pandas DataFrame containing the data
        :return: None
        """

        if not isinstance(data, DataFrame):
            raise TypeError(
                "There was an error reading the data: "
                "DataFrame has not been generated"
            )

        for x in ["time", "X", "experiments"]:
            if x not in data.columns:
                raise ValueError(f"Column {x} is missing from the dataset")
        if data.columns[0] != "experiments":
            raise ValueError("First column should be 'experiments'")
        if data.columns[1] != "time":
            raise ValueError("Second column should be 'time'")

        if len(data.columns) <= 3:
            raise ValueError(f"Data does not contain any metabolite columns")

        for x in data.columns:
            if x != "experiments" and data[x].dtypes != np.int64 and data[
                    x].dtypes != np.float64:
                raise ValueError(
                    f"Column {x} has values that are not of numeric type"
                )
            if all(data[x].isnull()) or all(data[x].isna()):
                raise ValueError(
                    f"The column {x} contains only null or NA values"
                )

        # To avoid errors when concatenating dataframes for the final summary
        data["experiments"] = data["experiments"].str.replace(pat=" ",
                                                              repl="_")

    @staticmethod
    def get_model_list():
        model_dir = Path(__file__).parent / "models"
        # Make fake dataframe
        df = DataFrame(
            columns=["time", "X"],
        )
        for file in os.listdir(str(model_dir)):
            if "model_" in file:
                module = importlib.import_module(
                    f"physiofit.models.{file[:-3]}"
                )
                model_class = getattr(module, "ChildModel")
                model = model_class(df)
                print(model.name)
        return

    def get_models(
            self,
            data: pd.DataFrame = None
    ):
        """
        Read modules containing the different models and add them to models
        attribute

        :return: list containing the different model objects
        """

        model_dir = Path(__file__).parent / "models"
        for file in os.listdir(str(model_dir)):
            if "model_" in file:
                module = importlib.import_module(
                    f"physiofit.models.{file[:-3]}"
                )
                model_class = getattr(module, "ChildModel")
                if data is not None:
                    self.models.append(model_class(data))
                else:
                    self.models.append(model_class(self.data))

    @staticmethod
    def get_local_model_folder() -> str:
        """
        Return the path towards the actual environment's used models folder
        """

        model_dir = Path(__file__).parent / "models"
        return str(model_dir)

    # TODO: Implement this function to add model to model folder and add
    #  button to GUI
    @staticmethod
    def add_model(model_file):
        pass

    @staticmethod
    def read_yaml(yaml_file: str | bytes) -> ConfigParser:

        """
        Import yaml configuration file and parse keyword arguments

        :param yaml_file: path to the yaml file or json file :return
        config_parser: Dictionary containing arguments parsed from yaml file

        """

        # Load config file
        try:
            if isinstance(yaml_file, str) or issubclass(type(yaml_file),
                                                        BytesIO):
                config_parser = ConfigParser.from_file(yaml_file)
            else:
                raise TypeError(
                    f"Trying to read object that is not a file or path to "
                    f"file: {yaml_file}"
                )
        except Exception as e:
            raise IOError(
                f"Error while reading yaml configuration file. "
                f"Traceback:{e}"
            )
        return config_parser

    def initialize_fitter(self, data: pd.DataFrame, **kwargs) -> PhysioFitter:
        """
        Initialize a PhysioFitter object

        :param data: input data
        :param kwargs: Keyword arguments for fitter initialization
        :return: None
        """
        #
        # wrong_keys = []

        # Initialize fitter
        fitter = PhysioFitter(
            data=data,
            model=kwargs["model"],
            mc=kwargs["mc"] if "mc" in kwargs else True,
            iterations=kwargs["iterations"] if "iterations" in kwargs else 100,
            sd=kwargs["sd"] if "sd" in kwargs else StandardDevs(),
            debug_mode=kwargs[
                "debug_mode"] if "debug_mode" in kwargs else False
        )

        if "sd" not in kwargs:
            fitter.sd.update(
                {"X": 0.2}
            )
            for col in self.data.columns[2:]:
                fitter.sd.update(
                    {col: 0.2}
                )

        fitter.initialize_sd_matrix()
        fitter.verify_attrs()

        return fitter

    def output_pdf(self, fitter: PhysioFitter, export_path: str | Path = None):
        """
        Handle the creation and output of a pdf file containing fit results as
        a plot

        :param fitter:
        :param export_path: Path to exported pdf. In local mode, it is sent to
                            the _res directory
        :return: None
        """

        if not self.figures:
            self.plot_data(fitter)

        try:
            with PdfPages(rf"{export_path}\plots.pdf") as pdf:
                for fig in self.figures:
                    pdf.savefig(fig[1])
        except Exception as e:
            raise IOError(f"Error while generating pdf:\n{e}")

    def output_plots(self, fitter, export_path):
        """
        Handle the creation and export of the different plots in svg format
        :return: None
        """

        if not self.figures:
            self.plot_data(fitter)

        try:
            for fig in self.figures:
                fig[1].savefig(rf"{export_path}\{fig[0]}.svg")
        except Exception:
            raise RuntimeError("Unknown error while generating output")

    def output_recap(self, export_path: str, galaxy=False):

        if not isinstance(self.multiple_experiments, list):
            raise TypeError(
                "The multiple experiments attribute must be a list"
            )
        if not self.multiple_experiments:
            raise ValueError(
                f"It seems that the multiple experiments list is empty: "
                f"{self.multiple_experiments}"
            )
        for idx, element in enumerate(self.multiple_experiments):
            if not isinstance(element, DataFrame):
                raise TypeError(
                    f"All the elements of multiple_experiments must be "
                    f"DataFrames. Wrong element type"
                    f"detected at indice {idx}"
                )
        final_df = concat(self.multiple_experiments)
        final_df = final_df.reset_index()
        final_df[["experiments", "parameter name"]] = final_df[
            "index"].str.split(" ", expand=True)
        final_df.set_index(["experiments", "parameter name"], inplace=True)
        final_df.drop("index", axis=1, inplace=True)
        if galaxy:
            final_df.to_csv(str(Path(export_path)))
        else:
            final_df.to_csv(f"{str(Path(export_path))}/summary.csv")

    @staticmethod
    def output_report(fitter, export_path: str | list = None):
        """
        Handle creation and export of the report containing stats from monte
        carlo analysis of optimization parameters

        :param export_path: Path to export the report
        :param fitter: PhysioFitter object containing results from the
                        optimization of parameters
        :param export_path: list of paths to export the stats and fluxes. [
                            0] is for stats and [1] for fluxes.
        """

        if isinstance(export_path, list):
            if len(export_path) != 2:
                raise ValueError(
                    f"Expected only 2 export paths and {len(export_path)}"
                    f" were detected"
                )
            stat_path = export_path[0]
            flux_path = export_path[1]
        else:
            flux_path = fr"{export_path}\flux_results.tsv"
            stat_path = fr"{export_path}\stat_results.tsv"

        # Get optimization results as dataframe
        opt_data = DataFrame.from_dict(fitter.parameter_stats,
                                       orient="columns")

        # Use IDs to clarify which parameter is described on each line
        opt_data.index = [param for param in fitter.model.parameters.keys()]
        opt_data.to_csv(flux_path, sep="\t")

        # Get khi² test results as dataframe and write to file
        with open(stat_path, "w+") as stat_out:
            if isinstance(fitter.khi2_res, DataFrame):
                stat_out.write("==================\n"
                               "Khi² test results\n"
                               "==================\n\n")
                stat_out.write(
                    fitter.khi2_res.to_string(
                        header=False, justify="center"
                    )
                )
                if fitter.khi2_res.at["p_val", "Values"] < 0.95:
                    stat_out.write(
                        f"\n\nAt level of 95% confidence, the model fits the "
                        f"data good enough with respect to the provided "
                        f"measurement SD. Value: "
                        f"{fitter.khi2_res.at['p_val', 'Values']}"
                    )

                else:
                    stat_out.write(
                        f"\n\nAt level of 95% confidence, the model does not "
                        f"fit the data good enough with respect to the "
                        f"provided measurement SD. Value: "
                        f"{fitter.khi2_res.at['p_val', 'Values']}\n"
                    )
            else:
                stat_out.write(
                    "No khi² test results available. The model has not been "
                    "tested against the data"
                )

            # Get AIC results as dataframe and write to file
            if isinstance(fitter.aic_res, DataFrame):
                stat_out.write(
                    "\n\n==================\n"
                    "AIC test results\n"
                    "==================\n\n"
                )
                stat_out.write(
                    fitter.aic_res.to_string(
                        header=False, justify="center"
                    )
                )
            else:
                stat_out.write(
                    "No AIC test results available."
                )

    def _get_plot_data(self, fitter):
        """
        Prepare data for plotting
        """

        # Initialize vectors and data for plotting
        if fitter.large_time_vector is not None:
            x = fitter.large_time_vector
        else:
            raise ValueError(
                "PhysioFitter model time_vector has not been initialized. "
                "Have you loaded in the data correctly?"
            )
        if fitter.experimental_matrix is not None:
            exp_mat = fitter.experimental_matrix
        else:
            raise ValueError(
                "PhysioFitter object does not have experimental data loaded in"
            )
        if fitter.large_matrix is not None:
            sim_mat = fitter.large_matrix
        else:
            raise ValueError("PhysioFitter simulated data does not exist yet")
        if fitter.matrices_ci is not None:
            sim_mat_ci = fitter.matrices_ci
            self.simulated_data_low_ci = DataFrame(
                columns=fitter.model.name_vector,
                index=x,
                data=sim_mat_ci["lower_ci"]
            )
            self.simulated_data_high_ci = DataFrame(
                columns=fitter.model.name_vector,
                index=x,
                data=sim_mat_ci["higher_ci"]
            )
        else:
            logger.warning(
                "Monte Carlo analysis has not been done, "
                "confidence intervals will not be computed"
            )

        self.experimental_data = DataFrame(
            columns=fitter.model.name_vector,
            index=fitter.model.time_vector,
            data=exp_mat
        ).sort_index(level=0)
        self.simulated_data = DataFrame(
            columns=fitter.model.name_vector,
            index=x,
            data=sim_mat
        ).sort_index(level=0)

    def plot_data(self, fitter, display: bool = False):
        """
        Plot the data

        :param display: Should plots be displayed or not on creation
        :param fitter: PhysioFitter object after optimization of parameters
        has been executed
        """

        self._get_plot_data(fitter)
        self._draw_plots(fitter, display)

    def _draw_plots(self, fitter, display: bool):
        """
        Draw plots and assign them to the _figures attribute for later access
        in pdf and plot generation functions

        :param display: Should plots be displayed or not on creation
        """

        for element in fitter.model.name_vector:

            # Initialize figure object
            fig, ax = plt.subplots()
            fig.set_size_inches(9, 5)

            # Draw experimental points onto the Axe
            exp = ax.scatter(
                self.experimental_data.index,
                self.experimental_data[element],
                marker='o',
                color="orange"
            )
            exp.set_label(f"Exp {element}")

            # Draw the simulated line onto the Axe
            sim_line, = ax.plot(
                self.simulated_data.index,
                self.simulated_data[element],
                linestyle='--'
            )
            sim_line.set_label(f"Sim {element}")

            # If monte carlo analysis has been done add the
            # corresponding sd to the plot as a red area
            if fitter.matrices_ci is not None:
                ax = self._add_sd_area(element, ax)

            # Finishing touches
            ax.set(xlim=0, ylim=0)
            ax.set_xlabel("Time", fontsize=21)
            ax.set_ylabel("Concentration", fontsize=21)
            ax.legend(prop={"size": 18})
            ax.set_title(f"{element}", fontsize=23)
            ax.tick_params(axis='both', which='major', labelsize=18)
            fig.tight_layout()

            # Add the figure with the metabolite name as index
            # to the _figures attribute for later use
            self.figures.append((element, fig))

        if display:
            for _ in self.figures:
                plt.show()
        plt.close()

    def _add_sd_area(self, element: str, ax: plt.Axes):
        """
        Draw red area around the fitting line to show confidence intervals

        :param element: Which variable is being plotted
        :param ax: axes on which to draw the area before
                   returning to mother function
        """

        y1 = self.simulated_data_low_ci[element].to_numpy()
        y2 = self.simulated_data_high_ci[element].to_numpy()
        x = self.simulated_data.index
        ax.fill_between(x, y1, y2, alpha=.3, linewidth=0, color="red")
        return ax


class ConfigParser:
    """
    The ConfigParser class is used to parse configuration files for the
    PhysioFit package. It reads a YAML file and extracts the necessary
    parameters for the model fitting process. It also includes methods to
    update the model with the parsed parameters and export the run
    parameters back to a yaml config file.
    """

    AUTHORIZED_KEYS = ["path_to_data", "model", "sds", "mc",
                       "iterations"]

    def __init__(
            self,
            selected_model,
            sds,
            mc,
            iterations,
            path_to_data=None
    ):

        self.path_to_data: str = path_to_data
        self.model: physiofit.models.base_model.Model = selected_model
        self.sds: physiofit.models.base_model.StandardDevs = StandardDevs(sds)
        self.mc: bool = mc
        self.iterations: int = iterations

        if not isinstance(self.mc, bool):
            raise TypeError(
                f"The MonteCarlo option must be given as a boolean (True or "
                f"False). Detected input: {self.mc},"
                f"type: {type(self.mc)}"
            )
        if not isinstance(self.iterations, int):
            raise TypeError(
                f"Number of iterations must be an integer: Detected input: "
                f"{self.mc}, type: {type(self.iterations)}"
            )

    @classmethod
    def from_file(cls, yaml_file):
        data = yaml.safe_load(yaml_file) if not isinstance(yaml_file, str) \
            else yaml.safe_load(open(yaml_file, 'r'))
        missing_keys = [key for key in cls.AUTHORIZED_KEYS if
                        key not in data.keys()]
        if missing_keys:
            raise ValueError(
                f"The keys {missing_keys} are missing from the input config "
                f"file")
        return cls(
            path_to_data=data["path_to_data"],
            selected_model=data["model"],
            sds=data["sds"],
            mc=data["mc"],
            iterations=data["iterations"]
        )

    @classmethod
    def from_galaxy(cls, galaxy_yaml):
        pass

    def get_kwargs(self):
        return {
            "path_to_data": self.path_to_data,
            "model": self.model,
            "mc": self.mc,
            "iterations": self.iterations,
            "sd": self.sds
        }

    def check_data_path(self):
        """
        Check if the data path is valid
        :return: True if the path is valid, False otherwise, None if no path
        """
        if self.path_to_data:
            return Path(self.path_to_data).is_file()
        else:
            raise ValueError("No data path has been provided")

    def update_model(self, model):

        if self.model["parameters_to_estimate"] is not None:
            model.parameters.update(self.model["parameters_to_estimate"])
        if self.model["bounds"] is not None:
            model.bounds.update(Bounds(self.model["bounds"]))
        return model

    def export_config(self, export_path):

        with open(fr"{export_path}/config_file.yml", "w") as file:
            data = {
                "model": {
                    "model_name": self.model.name,
                    "parameters_to_estimate": self.model.parameters,
                    "bounds": {name: f"{bounds[0], bounds[1]}" for name, bounds
                               in self.model.bounds.items()},
                    "args": self.model.args,
                },
                "sds": dict(self.sds),
                "mc": self.mc,
                "iterations": self.iterations,
                "path_to_data": str(self.path_to_data)
            }
            yaml.safe_dump(
                data,
                file
            )
