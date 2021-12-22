"""Module to handle inputs and outputs for the Physiofit software"""
from pathlib import Path

from pandas import DataFrame, read_csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from physiofit.base.fitter import PhysioFitter


class IoHandler:

    allowed_keys = {"vini", "mc", "iterations", "pos", "conc_biom_bounds", "flux_biom_bounds", "conc_met_bounds",
                    "flux_met_bounds", "weight", "sd_X", "sd_M", "save", "summary", "debug_mode"}

    def __init__(self, source='local'):
        """
        Input/Output class that handles the former and initializes the PhysioFitter component object. It is the
        interface for interacting with the Physiofit package.
        Usage:  iohandle = IoHandler(source)
                iohandle.local_in(data, kwargs) --> kwargs are passed on to the PhysioFitter class for initialization
                iohandle.fitter.fitter_function()
                ...
                iohandle.local_out(*args) --> args define the type of output that should be generated

        :param source: Input source
        """

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

    def local_in(self, data, **kwargs):
        """Function for reading data and initializing the fitter object in local context"""

        if self.data:
            raise ValueError(f"It seems data is already loaded in. Data: \n{self.data}\nHome path: {self.home_path}")
        elif self.input_source == "local":
            self.home_path = Path(data).parent.resolve()
            if not self.home_path.exists():
                raise KeyError(f"The input file does not exist. Path: {self.home_path}")
            self.data = IoHandler._read_data(data)
            self.data = self.data.sort_values("time", ignore_index=True)
            self.names = self.data.columns[1:].to_list()
        else:
            raise IOError(f"Wrong input source selected. Source: {self.input_source}")
        self.initialize_fitter(kwargs)

    def local_out(self, *args):
        """Function for coordinating exports in local context"""

        for arg in args:
            if arg not in ["data", "plot", "pdf"]:
                raise ValueError(f"Detected argument is not an output type: {arg}.\n"
                                 f"Accepted output types are: data, plot, pdf")
            else:
                self._output_type.append(arg)
        if "data" in self._output_type:
            self._output_report()
        if "plot" in self._output_type:
            self._output_plots()
        if "pdf" in self._output_type:
            self._output_pdf()

    def initialize_fitter(self, kwargs):
        """
        Initialize the PhysioFitter object and pass kwargs on down
        :param kwargs: Keyword arguments for fitter initialization
        :return: None
        """

        wrong_keys = []
        self.fitter = PhysioFitter(self.data)
        self.fitter.__dict__.update((key, value) if key in IoHandler.allowed_keys else wrong_keys.append(key)
                                    for key, value in kwargs.items())
        # TODO: Think of how to implement this: should we use properties instead?
        self.fitter.initialize_vectors()
        self.fitter.initialize_weight_matrix()
        self.fitter.initialize_bounds()
        self.fitter.logger.debug(f"Fitter attribute dictionnary:\n{self.fitter.__dict__}")
        if wrong_keys:
            raise KeyError(f"Some keyword arguments were not valid: {wrong_keys}")

    def _output_pdf(self):
        """Handle the creation and output of a pdf file containing fit results in plot form"""

        if not self.home_path:
            raise RuntimeError("No home path detected. Was data loaded in correctly?")
        if not self._figures:
            raise RuntimeError("Plots have not been created. Please launch the plot_data() function first")
        try:
            with PdfPages(rf"{self.home_path}\plots.pdf") as pdf:
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
                fig[1].savefig(rf"{self.home_path}\{fig[0]}.svg")
        except Exception:
            raise RuntimeError("Unknown error while generating output")

    def _output_report(self):
        """
        Handle creation and export of the report containing stats from monte carlo analysis of optimization
        parameters
        """

        self.fitter.logger.debug(f"Parameter Stats:\n{self.fitter.parameter_stats}")
        opt_data = DataFrame.from_dict(self.fitter.parameter_stats,
                                       orient="columns")
        # Use IDs to clarify which parameter is described on each line
        opt_data.index = self.fitter.ids
        opt_data.to_csv(fr"{self.home_path}\Optimized_parameter_statistics.csv")

    def _get_plot_data(self):
        """
        Prepare data for plotting
        """

        if hasattr(self.fitter, "time_vector"):
            x = self.fitter.time_vector
        else:
            raise ValueError("PhysioFitter time_vector has not been initialized. "
                             "Have you loaded in the data correctly?")
        if hasattr(self.fitter, "experimental_matrix"):
            exp_mat = self.fitter.experimental_matrix
        else:
            raise ValueError("PhysioFitter object does not have experimental data loaded in")
        if hasattr(self.fitter, "simulated_matrix"):
            sim_mat = self.fitter.simulated_matrix
        else:
            raise ValueError("PhysioFitter simulated data does not exist yet")
        if hasattr(self.fitter, "matrices_sds"):
            sim_mat_sds = self.fitter.matrices_sds
        else:
            raise ValueError("PhysioFitter monte carlo analysis has not been done")
        self.experimental_data = DataFrame(columns=self.names, index=x, data=exp_mat)
        self.simulated_data = DataFrame(columns=self.names, index=x, data=sim_mat)
        self.simulated_data_sds = DataFrame(columns=self.names, index=x, data=sim_mat_sds)

    def plot_data(self, display=True):
        """
        Plot the extracellular flux data
        :param display: should plots be displayed
        """

        self._get_plot_data()
        self._draw_plots(display)

    def _draw_plots(self, display):
        """
        Draw the plots and assign them to the _figures attribute for later acces in pdf and plot generation functions

        :param display: Should plots be displayed or not on creation
        """

        for element in self.names:
            fig, ax = plt.subplots()
            fig.set_size_inches(18.5, 10.5)
            exp = ax.scatter(self.experimental_data.index,
                             self.experimental_data[element],
                             marker='o', color="orange")
            exp.set_label(f"Exp {element}")
            sim_line, = ax.plot(self.simulated_data.index,
                                self.simulated_data[element],
                                linestyle='--')
            sim_line.set_label(f"Sim {element}")
            ax = self._add_sd_area(element, ax)
            ax.set(xlim=0, ylim=0, xlabel="Time", ylabel="Concentration")
            ax.legend()
            ax.set_title(f"{element} extracellular flux")
            fig.tight_layout()
            self._figures.append((element, fig))
        if display:
            for _ in self._figures:
                plt.show()

    def _add_sd_area(self, element, ax):
        """
        Draw red area around the fitting line to show SD

        :param element: Which variable is being plotted
        :param ax: axes on which to draw the area before returning to mother function
        """

        y1 = self.simulated_data[element].to_numpy() + self.simulated_data_sds[element].to_numpy()
        y2 = self.simulated_data[element].to_numpy() - self.simulated_data_sds[element].to_numpy()
        x = self.simulated_data.index
        ax.fill_between(x, y1, y2, alpha=.3, linewidth=0, color="red")
        return ax

    @staticmethod
    def _read_data(path_to_data: str) -> DataFrame:
        """
        Read initial data file (csv or tsv)

        :param path_to_data: str containing the relative or absolute path to the data
        :return: pandas DataFrame containing the data
        """

        data_path = Path(path_to_data).resolve()
        if data_path.suffix == ".tsv":
            data = read_csv(data_path, sep="\t")
        elif data_path.suffix == ".csv":
            data = read_csv(data_path, sep=";")
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
                raise ValueError(f"The column {col} is missing from the dataset")
        if len(data.columns) <= 2:
            raise ValueError(f"The data does not contain any metabolite columns")
        for col in data.columns:
            if data[col].dtypes != np.int64 and data[col].dtypes != np.float64:
                raise ValueError(f"The column {col} has values that are not of numeric type")


if __name__ == "__main__":
    iostream = IoHandler()
    iostream.local_in(
        r"C:\Users\legregam\Documents\Projets\PhysioFit\Example\KEIO_test_data\KEIO_ROBOT6_7.tsv",
        iterations=50, vini=0.05, weight=[0.02, 0.46, 0.1]
    )
    iostream.fitter.optimize()
    iostream.fitter.monte_carlo_analysis()
    iostream.plot_data(False)
    iostream.local_out("data", "plot", "pdf")
