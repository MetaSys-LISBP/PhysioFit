"""Module to handle inputs and outputs for the Physiofit software"""
from pathlib import Path

from pandas import DataFrame, read_csv
import numpy as np
import matplotlib.pyplot as plt

from physiofit.base.fitter import PhysioFitter


class IoHandler:

    def __init__(self, source='local'):

        self._input_source = source
        self._output_type = "data"
        self.output = None
        self.data = None
        self.names = None
        self.simulated_data_sds = None
        self.simulated_data = None
        self.experimental_data = None
        self.home_path = None

    def handle_local_data(self, input_):

        if self._input_source == "local":
            if self._output_type == "data":
                self.home_path = Path(input_).parent.resolve()
                self.data = IoHandler._read_data(input_)
                self.data = self.data.sort_values("time", ignore_index=True)
                self.names = self.data.columns[1:].to_list()
            elif self._output_type == "plot":
                if not self.home_path:
                    raise RuntimeError("No home path detected. Was data loaded in correctly?")
                try:
                    res_dir = self.home_path / "Results"
                    res_dir.mkdir()
                    for fig in input_:
                        fig[1].savefig(rf"{res_dir}\{fig[0]}.svg")
                except TypeError:
                    raise
                except FileExistsError:
                    raise("The Results folder already exists. Have you already performed this run?")
                except Exception:
                    raise RuntimeError("Unknown error while generating output")
        else:
            raise IOError(f"Wrong input source selected. Source: {self._input_source}")

    def _get_plot_data(self, fitter):

        x = fitter.time_vector
        exp_mat = fitter.experimental_matrix
        sim_mat = fitter.simulated_matrix
        self.experimental_data = DataFrame(columns=self.names, index=x, data=exp_mat)
        self.simulated_data = DataFrame(columns=self.names, index=x, data=sim_mat)
        self.simulated_data_sds = fitter.matrices_sds

    def _plot_pdf(self):
        pass

    def plot_data(self, fitter, pdf=False, display=True):

        self._get_plot_data(fitter)
        if pdf:
            self._plot_pdf()
        else:
            self._draw_plots(display)

    def _draw_plots(self, display):

        figures = []
        self._output_type = "plot"
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
            ax.set(xlim=0, ylim=0, xlabel="Time", ylabel="Concentration")
            ax.legend()
            ax.set_title(f"{element} extracellular flux")
            fig.tight_layout()
            figures.append((element, fig))
        if display:
            for _ in figures:
                plt.show()
        if self._input_source == "local":
            self.handle_local_data(figures)

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
    iostream.handle_local_data(
        r"C:\Users\legregam\Documents\Projets\PhysioFit\Example\KEIO_test_data\KEIO_ROBOT6_7.tsv"
    )
    test = PhysioFitter(iostream.data, vini=0.05, weight=[0.02, 0.46, 0.1])
    test.optimize()
    test.monte_carlo_analysis()
    iostream.plot_data(test)
