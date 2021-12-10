"""
PhysioFit software main module
"""

from pathlib import Path
import logging

import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame
from scipy.optimize import minimize

from physiofit.logger import initialize_fitter_logger

mod_logger = logging.getLogger("PhysioFit.base.fitter")


class PhysioFitter:

    def __init__(self, path_to_data, vini=0.04, mc=True, iterations=50, pos=True, conc_biom_bounds=(1e-2, 50),
                 flux_biom_bounds=(0.01, 50), conc_met_bounds=(1e-6, 50), flux_met_bounds=(-50, 50), weight=None,
                 sd_X=0.002, sd_M=0.5, save=True, summary=False, debug_mode=True):

        """
        The PhysioFitter class is responsible for most of Physiofit's heavy lifting. Features included are:
            * loading of data from csv or tsv file
            * equation system initialization using the following analytical functions (in absence of lag and
              degradation:
                X(t) = X0 * exp(mu * t)
                Mi(t) = qMi * (X0 / mu) * (exp(mu * t) - 1) + Mi0
            * simulation of data points from given initial parameters
            * cost calculation using the equation:
                residuum = sum((sim - meas) / weight)²
            * optimization of the initial parameters using scipy.optimize.minimize ('L-BFGS-B' method)
            * calling the Stat Analyzer objet for sensitivity analysis, khi2 tests and plotting (see documentation
              relative to the component Stat Analyzer class for more details)

        :param path_to_data: path to input file
        :param mc: Should Monte-Carlo sensitivity analysis be performed (default=True)
        :param vini: initial value for fluxes and concentrations (default=1)
        :type vini: int or float
        :type mc: Boolean
        :param iterations: number of iterations for Monte-Carlo simulation (default=50)
        :type iterations: int
        :param pos: Negative concentrations of noisy datasets generated during Monte-Carlo iterations are set to 0 if
                    True (default=True)
        :type pos: Boolean
        :param conc_biom_bounds: lower and upper bounds for biomass concentration (X0)
        :type conc_biom_bounds: tuple of ints/floats (lower, upper)
        :param flux_biom_bounds: lower and upper bounds for biomass rate of change (mu)
        :type flux_biom_bounds: tuple of ints/floats (lower, upper)
        :param conc_met_bounds: lower and upper bounds for metabolite concentration (mi0)
        :type conc_met_bounds: tuple of ints/floats (lower, upper)
        :param flux_met_bounds: lower and upper bounds for metabolite rate of change (qi)
        :type flux_met_bounds: tuple of ints/floats (lower, upper)
        :param weight: weight matrix used for residuum calculations. Can be:
                       * a matrix with the same dimensions as the measurements matrix (but without the time column)
                       * a named vector containing weights for all the metabolites provided in the input file
                       * 0 (by default), in which case the matrix is automatically loaded from the file xxx_sd.csv/.tsv
                       (where xxx is the data file name) if the file exists. Otherwise, weight is constructed from sd_X
                       and sd_M arguments
        :type weight: int, float, list or ndarray
        :param sd_X: Standard deviation on biomass concentrations (default = 0.002), used only if weight = 0
        :type sd_X: int or flaot
        :param sd_M: Standard deviation on metabolite concentrations (defaul = 0.5), used only if weight = 0
        :type sd_M: int or float
        :param save: Should results be saved
        :type save: Boolean
        :param summary: Should results of the khi-2 test be displayed
        :type summary: Boolean
        """

        self.data = PhysioFitter._read_data(path_to_data)
        self.data = self.data.sort_values("time", ignore_index=True)
        self.vini = vini
        self.mc = mc
        self.iterations = iterations
        self.pos = pos
        self.conc_biom_bounds = conc_biom_bounds
        self.flux_biom_bounds = flux_biom_bounds
        self.conc_met_bounds = conc_met_bounds
        self.flux_met_bounds = flux_met_bounds
        self.weight = weight
        self.sd_X = sd_X
        self.sd_M = sd_M
        self.save = save
        self.summary = summary
        self.debug_mode = debug_mode
        self.logger = initialize_fitter_logger(self.debug_mode)

        self.simulated_matrix = None
        self.optimize_results = None
        self.params = None
        self.ids = None
        self.bounds = None

        self.logger.info("Initializing vectors...\n")
        self._initialize_vectors()
        self.logger.debug(f"Time vector: {self.time_vector}\n"
                          f"Name vector: {self.name_vector}\n"
                          f"Experimental Data: \n{self.data}\n"
                          f"Parameters: {self.ids}\n"
                          f"Parameter vector: {self.params}\n")
        self.logger.info("Initializing Weight matrix...\n")
        self._initialize_weight_matrix()
        self.logger.debug(f"Weight matrix: {self.weight}")
        self.logger.info("Initializing bounds...")
        self._initialize_bounds()
        self.logger.debug(f"Bounds: {self.bounds}")

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
        PhysioFitter._verify_data(data)
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

    def _initialize_vectors(self):
        """
        Initialize the vectors needed for flux calculations from the input parameters

        :return: None
        """

        self.time_vector = self.data.time.to_numpy()
        self.name_vector = self.data.drop("time", axis=1).columns.to_list()
        self.exp_data_matrix = self.data.drop("time", axis=1).to_numpy()
        metabolites = self.name_vector[1:]
        mu = self.vini
        x_0 = self.vini
        self.params = [x_0, mu]
        self.ids = ["X_0", "mu"]
        for met in metabolites:
            self.params.append(self.vini)
            self.params.append(self.vini)
            self.ids.append(f"{met}_q")
            self.ids.append(f"{met}_M0")

    def _initialize_weight_matrix(self):
        """
        Initialize the weight matrix from different types of inputs: single value, vector or matrix.

        :return: None
        """

        # TODO: This function can be optimized, if the input is a matrix we should detect it directly

        # When 0 is given as input weight, we assume the weights are given in an external file
        if self.weight is None:
            self._read_weight_file()
            return
        # When weight is a single value, we build a weight matrix containing the value in all positions
        if isinstance(self.weight, int) or isinstance(self.weight, float):
            self._build_weight_matrix()
            return
        if not isinstance(self.weight, np.ndarray):
            if not isinstance(self.weight, list):
                raise TypeError(f"Cannot coerce weights to array. Please check that a list or array is given as input."
                                f"\nCurrent input: \n{self.weight}")
            else:
                self.weight = np.array(self.weight)
        if not np.issubdtype(self.weight.dtype, np.number):
            try:
                self.weight = self.weight.astype(float)
            except ValueError:
                raise ValueError(f"The weight vector/matrix contains values that are not numbers. \n"
                                 f"Current weight vector/matrix: \n{self.weight}")
            except Exception as e:
                raise RuntimeError(f"Unknown error: {e}")
        else:
            # If the array is not the right shape, we assume it is a vector that needs to be tiled into a matrix
            if self.weight.shape != self.exp_data_matrix.shape:
                try:
                    self._build_weight_matrix()
                except ValueError:
                    raise
                except RuntimeError:
                    raise
            else:
                return

    def _initialize_bounds(self):

        # We set the bounds for x0 and mu
        bounds = [
            self.conc_biom_bounds,
            self.flux_biom_bounds
        ]
        # We get the number of times that we must add the m0 and q0 bounds (once per metabolite)
        ids_range = int((len(self.ids) - 2) / 2)  # We force int so that Python does not think it could be float
        for _ in range(ids_range):
            bounds.append(
                self.flux_met_bounds  # q_0
            )
            bounds.append(
                self.conc_met_bounds  # M_0
            )
        self.bounds = tuple(bounds)

    def _read_weight_file(self):
        pass

    def _build_weight_matrix(self):
        """
        Build the weight matrix from different input types

        :return: None
        """

        if isinstance(self.weight, np.ndarray):
            if self.weight.size != self.exp_data_matrix[0].size:
                raise ValueError("Weight vector not of right size")
            else:
                self.weight = np.tile(self.weight, (len(self.exp_data_matrix), 1))
        elif isinstance(self.weight, int) or isinstance(self.weight, float):
            self.weight = np.full(self.exp_data_matrix.shape, self.weight)
        else:
            raise RuntimeError("Unknown error")

    def simulate(self, equation_type="simple"):

        if equation_type == "simple":
            self.simulated_matrix = PhysioFitter._simple_sim(self.params, self.exp_data_matrix, self.time_vector)
        else:
            pass

    def optimize(self):

        self.logger.info("\nRunning optimization...")
        self.optimize_results = PhysioFitter._run_optimization(self.params, self.exp_data_matrix,
                                                               self.time_vector, self.weight, self.bounds)
        self.logger.info(f"Optimization results: \n{self.optimize_results}")
        for i, param in zip(self.ids, self.optimize_results.x):
            self.logger.info(f"\n{i} = {param}")
        self.simulated_matrix = PhysioFitter._simple_sim(self.optimize_results.x, self.exp_data_matrix,
                                                         self.time_vector)
        self.logger.info(f"Final Simulated Matrix: \n{self.simulated_matrix}")

    def test_plot(self):

        fig, (ax1, ax2, ax3) = plt.subplots(3)
        fig.set_size_inches(18.5, 10.5)
        x = self.time_vector
        exp_biomass = self.exp_data_matrix[:, 0]
        exp_glc = self.exp_data_matrix[:, 1]
        exp_ace = self.exp_data_matrix[:, 2]
        sim_biomass = self.simulated_matrix[:, 0]
        sim_glc = self.simulated_matrix[:, 1]
        sim_ace = self.simulated_matrix[:, 2]

        # Work on biomass ax
        exp_biomass = ax1.scatter(x, exp_biomass, marker='o', color="orange")
        exp_biomass.set_label("Exp Biomass")
        sim_biomass_line, = ax1.plot(x, sim_biomass, linestyle='--')
        sim_biomass_line.set_label("Sim Biomass")
        ax1.set(xlim=0, ylim=0, xlabel="Time", ylabel="Concentration")
        ax1.legend()
        ax1.set_title("Sim/Exp comparison")

        # Work on glucose ax
        exp_glucose = ax2.scatter(x, exp_glc, marker='o', color="orange")
        exp_glucose.set_label("Exp Glucose")
        sim_glucose_line, = ax2.plot(x, sim_glc, linestyle='--')
        sim_glucose_line.set_label("Sim Glucose")
        ax2.set(xlim=0, ylim=0, xlabel="Time", ylabel="Concentration")
        ax2.legend()

        # Work on acetate ax
        exp_ace = ax3.scatter(x, exp_ace, marker='o', color="orange")
        exp_ace.set_label("Exp Acetate")
        sim_ace_line, = ax3.plot(x, sim_ace, linestyle='--')
        sim_ace_line.set_label("Sim Acetate")
        ax3.set(xlim=0, ylim=0, xlabel="Time", ylabel="Concentration")
        ax3.legend()

        plt.show()

    @staticmethod
    def _simple_sim(params, exp_data_matrix, time_vector):

        simulated_matrix = np.empty_like(exp_data_matrix)
        x_0 = params[0]
        mu = params[1]
        exp_mu_t = np.exp(mu * time_vector)
        simulated_matrix[:, 0] = x_0 * exp_mu_t
        for i in range(1, int(len(params) / 2)):
            q = params[i * 2]
            m_0 = params[i * 2 + 1]
            simulated_matrix[:, i] = q * (x_0 / mu) * (exp_mu_t - 1) + m_0
        return simulated_matrix

    @staticmethod
    def _calculate_cost(params, exp_data_matrix, time_vector, weight_matrix):

        simulated_matrix = PhysioFitter._simple_sim(params, exp_data_matrix, time_vector)
        cost_val = np.square((simulated_matrix - exp_data_matrix) / weight_matrix)
        residuum = np.nansum(cost_val)
        return residuum

    @staticmethod
    def _run_optimization(params, exp_data_matrix, time_vector, weight_matrix, bounds):

        optimize_results = minimize(PhysioFitter._calculate_cost, x0=params, args=(
            exp_data_matrix, time_vector, weight_matrix), method="L-BFGS-B", bounds=bounds)
        return optimize_results

    def monte_carlo_analysis(self):
        if not self.optimize_results:
            raise RuntimeError("Running Monte Carlo simulation without having run the optimization is impossible "
                               "as best fit results are needed to generate the simulated matrix")

        opt_res = self.optimize_results
        opt_params = []
        matrices = []
        for _ in range(self.iterations):
            new_matrix = self._apply_noise()
            opt_res = PhysioFitter._run_optimization(opt_res.x, new_matrix, self.time_vector,
                                                     self.weight, self.bounds)
            matrices.append(PhysioFitter._simple_sim(opt_res.x, new_matrix, self.time_vector))
            opt_params.append(opt_res.x)
        opt_params = np.array(opt_params)
        matrices = np.array(matrices)
        print(f"optimized parameters:\n {np.std(opt_params, 0)}")
        print(f"new matrices:\n {np.std(matrices, 0)}")
        return

    @staticmethod
    def _add_noise(array, sd):
        """
        Add random noise to a given array using input standard deviations

        :param array: input array on which to apply noise
        :type array: class: numpy.ndarray
        :param sd: standard deviation to apply to the input array
        :type sd: class: numpy.ndarray
        :return: noisy ndarray
        """

        output = np.random.normal(loc=array, scale=sd, size=array.size)
        return output

    def _apply_noise(self):
        """
        Apply noise to the simulated matrix obtained using optimized parameters. SDs are obtained from the weight
        matrix
        """

        new_matrix = np.array([
            PhysioFitter._add_noise(self.simulated_matrix[:, idx], sd)
            for idx, sd in enumerate(self.weight[0, :])
        ])

        return new_matrix.transpose()


if __name__ == "__main__":
    test = PhysioFitter(r"C:\Users\legregam\Documents\Projets\PhysioFit\Example\KEIO_test_data\KEIO_ROBOT6_7.tsv",
                        vini=0.05, weight=[0.02, 0.46, 0.1])
    test.optimize()
    # test.test_plot()
    print(test.simulated_matrix)
    print(test.optimize_results)
    test.monte_carlo_analysis()
    # out = PhysioFitter._add_noise(test.simulated_matrix[:, 1], 0.46)
    # print(out)
    # out = test.monte_carlo_analysis()
    # print(out)
    # print(out.shape)
    # print(out.ndim)
    # print(f"Simulated Matrix:\n {test.simulated_matrix}\n"
    #       f"STD results: \n {test.monte_carlo_analysis()}")
