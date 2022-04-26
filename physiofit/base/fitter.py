"""
PhysioFit software main module
"""

import logging

import numpy as np
from pandas import DataFrame
from scipy.optimize import minimize
from scipy.stats import chi2

from physiofit.logger import initialize_fitter_logger

mod_logger = logging.getLogger("PhysioFit.base.fitter")


# TODO: add estimate deg function (eq 6) with plot of best fit and measured values


class PhysioFitter:
    """
    This class is responsible for most of Physiofit's heavy lifting. Features included are:

        * loading of data from **csv** or **tsv** file
        * **equation system initialization** using the following analytical functions (in absence of lag and
          degradation:

            X(t) = X0 * exp(mu * t)
            Mi(t) = qMi * (X0 / mu) * (exp(mu * t) - 1) + Mi0

        * **simulation of data points** from given initial parameters
        * **cost calculation** using the equation:

            residuum = sum((sim - meas) / weight)²

        * **optimization of the initial parameters** using `scipy.optimize.minimize ('L-BFGS-B' method) <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html#optimize-minimize-lbfgsb>`_
        * **sensitivity analysis, khi2 tests and plotting**

    :param data: DataFrame containing data and passed by IOstream object
    :type data: class: pandas.DataFrame
    :param vini: initial value for fluxes and concentrations (default=1)
    :type vini: int or float
    :param mc: Should Monte-Carlo sensitivity analysis be performed (default=True)
    :type mc: Boolean
    :param iterations: number of iterations for Monte-Carlo simulation (default=50)
    :type iterations: int
    :param conc_biom_bounds: lower and upper bounds for biomass concentration (X0)
    :type conc_biom_bounds: tuple of ints/floats (lower, upper)
    :param flux_biom_bounds: lower and upper bounds for biomass rate of change (mu)
    :type flux_biom_bounds: tuple of ints/floats (lower, upper)
    :param conc_met_bounds: lower and upper bounds for metabolite concentration (mi0)
    :type conc_met_bounds: tuple of ints/floats (lower, upper)
    :param flux_met_bounds: lower and upper bounds for metabolite rate of change (qi)
    :type flux_met_bounds: tuple of ints/floats (lower, upper)
    :param sd: weight matrix used for residuum calculations. Can be:

                * a matrix with the same dimensions as the measurements matrix (but without the time column)
                * a named vector containing weights for all the metabolites provided in the input file
                * 0  in which case the matrix is automatically constructed from default values
                * a dictionary with the data column headers as keys and the associated value as a scalar or list

    :type sd: int, float, list, dict or ndarray
    :param deg: dictionary of degradation constants for each metabolite
    :type deg: dict
    :param t_lag: Should lag phase length be estimated
    :type t_lag: bool
    """

    def __init__(self, data, vini=0.2, mc=True, iterations=100, conc_biom_bounds=(1e-4, 50),
                 flux_biom_bounds=(1e-3, 2), conc_met_bounds=(1e-6, 50), flux_met_bounds=(-50, 50), sd=None,
                 deg=None, t_lag=False, debug_mode=False):

        self.data = data
        self.vini = vini
        self.mc = mc
        self.iterations = iterations
        self.conc_biom_bounds = conc_biom_bounds
        self.flux_biom_bounds = flux_biom_bounds
        self.conc_met_bounds = conc_met_bounds
        self.flux_met_bounds = flux_met_bounds
        self.sd = sd
        self.deg = deg
        self.t_lag = t_lag
        self.debug_mode = debug_mode
        if not hasattr(self, "logger"):
            self.logger = initialize_fitter_logger(self.debug_mode)

        self.simulated_matrix = None
        self.simulated_data = None
        self.optimize_results = None
        self.simulate = None
        self.time_vector = None
        self.name_vector = None
        self.deg_vector = None
        self.experimental_matrix = None
        self.params = None
        self.ids = None
        self.bounds = None
        self.parameter_stats = None
        self.opt_params_sds = None
        self.matrices_ci = None
        self.opt_conf_ints = None
        self.khi2_res = None

    def initialize_vectors(self):
        """
        Initialize the vectors needed for flux calculations from the input parameters

        :return: None
        """

        self.time_vector = self.data.time.to_numpy()
        self.name_vector = self.data.drop("time", axis=1).columns.to_list()
        self.experimental_matrix = self.data.drop("time", axis=1).to_numpy()
        metabolites = self.name_vector[1:]
        mu = self.vini
        x_0 = self.vini

        # Check if t_lag should be added to the list of parameters to estimate
        if self.t_lag:
            t_lag = 0.1 * self.time_vector.max()
            self.params = [x_0, mu, t_lag]
            self.ids = ["X_0", "mu", "t_lag"]
        else:
            self.params = [x_0, mu]
            self.ids = ["X_0", "mu"]

        # Handle the creation of a vector of degradation constants
        if self.deg:
            if type(self.deg) is dict:
                for key in self.deg.keys():
                    if key not in metabolites:
                        raise KeyError(f"The degradation constant for {key} is missing. If no degradation applies to this "
                                       f"metabolite, please enter 0 in the corresponding dictionary entry")
            self.deg_vector = [self.deg[met] for met in metabolites]
        elif self.deg == {} or self.deg is None:
            self.deg_vector = [0 for _ in metabolites]

        # Build a list containing each metabolite's q and m0
        for met in metabolites:
            self.params.append(self.vini)
            self.params.append(self.vini)
            self.ids.append(f"{met}_q")
            self.ids.append(f"{met}_M0")

    def verify_attrs(self):
        """Check that attributes are valid"""

        allowed_vinis = [int, float]
        if type(self.vini) not in allowed_vinis:
            raise TypeError(f"Initial value for fluxes and concentrations is not a number. Detected type:  "
                            f"{type(self.vini)}\nValid types: {allowed_vinis}")

        for bound in [self.conc_biom_bounds, self.flux_biom_bounds, self.conc_met_bounds, self.flux_met_bounds]:
            if type(bound) is not tuple:
                raise TypeError(f"Error reading bounds. Bounds should be ('int/float', 'int/float') tuples.\n"
                                f"Current bounds: {bound}")
            if self.vini < bound[0] or self.vini > bound[1]:
                raise RuntimeError(f"Initial value for fluxes and concentrations ({self.vini}) cannot be set outside "
                                   f"the given bounds: {bound}")

        if type(self.iterations) is not int:
            raise TypeError(f"Number of monte carlo iterations must be an integer, and not of type "
                            f"{type(self.iterations)}")

        allowed_weights = [int, float, list, np.ndarray]
        if type(self.sd) not in allowed_weights:
            raise TypeError(f"Weights is not in the right format ({type(self.sd)}. "
                            f"Compatible formats are:\n{allowed_weights}")

        if type(self.deg_vector) is not list:
            raise TypeError(f"Degradation constants have not been well initialized.\nConstants: {self.deg}")

        if type(self.t_lag) is not bool:
            raise TypeError(f"t_lag parameter must be a boolean (True or False)")

    def _weight_dict_to_matrix(self):
        """Convert weight dictionary to matrix/vector"""

        # Perform checks
        for key in self.sd.keys():
            if key not in self.name_vector:
                raise KeyError(f"The key {key} is not part of the data headers")
        for name in self.name_vector:
            if name not in self.sd.keys():
                raise KeyError(f"The key {name} is missing from the weights dict")

        # Get lengths of each weight entry
        weight_lengths = [len(self.sd[key]) if type(self.sd[key]) not in [float, int] else 1
                          for key in self.sd.keys()]

        # Make sure that lengths are the same
        if not all(elem == weight_lengths[0] for elem in weight_lengths):
            raise ValueError("All weight vectors must have the same length")

        # Build matrix/vector
        if weight_lengths[0] == 1:
            self.sd = [self.sd[name] for name in self.name_vector]
        else:
            columns = (self.sd[name] for name in self.name_vector)
            matrix = np.column_stack(columns)
            self.sd = matrix

    def initialize_weight_matrix(self):
        """
        Initialize the weight matrix from different types of inputs: single value, vector or matrix.

        :return: None
        """

        # This function can be optimized, if the input is a matrix we should detect it directly
        self.logger.info("Initializing Weight matrix...\n")

        # If weight is None, we generate the default matrix
        if self.sd is None:
            try:
                self.sd = {"X": 0.02}
                for col in self.data.columns[2:]:
                    self.sd.update({col: 0.5})
            except Exception:
                raise

        if type(self.sd) is dict:
            self._weight_dict_to_matrix()
        # When weight is a single value, we build a weight matrix containing the value in all positions
        if isinstance(self.sd, int) or isinstance(self.sd, float):
            self._build_weight_matrix()
            self.logger.debug(f"Weight matrix: {self.sd}\n")
            return
        if not isinstance(self.sd, np.ndarray) and not isinstance(self.sd, list):
            raise TypeError(
                f"Cannot coerce weights to array. Please check that a list or array is given as input."
                f"\nCurrent input: \n{self.sd}")
        else:
            self.sd = np.array(self.sd)
        if not np.issubdtype(self.sd.dtype, np.number):
            try:
                self.sd = self.sd.astype(float)
            except ValueError:
                raise ValueError(f"The weight vector/matrix contains values that are not numeric. \n"
                                 f"Current weight vector/matrix: \n{self.sd}")
            except Exception as e:
                raise RuntimeError(f"Unknown error: {e}")
        else:
            # If the array is not the right shape, we assume it is a vector that needs to be tiled into a matrix
            if self.sd.shape != self.experimental_matrix.shape:
                try:
                    self._build_weight_matrix()
                except ValueError:
                    raise
                except RuntimeError:
                    raise
            else:
                self.logger.debug(f"Weight matrix: {self.sd}\n")
                return
        self.logger.info(f"Weight Matrix:\n{self.sd}\n")

    def initialize_equation(self):

        if self.t_lag and self.deg:
            self.logger.debug("_total_sim function used for simulation\n")
            self.simulate = self._total_sim
        if self.t_lag and not self.deg:
            self.logger.debug("_lag_sim function used for simulation\n")
            self.simulate = self._lag_sim
        if not self.t_lag:
            self.logger.debug("_simple_sim function used for simulation\n")
            self.simulate = self._simple_sim

    def initialize_bounds(self):
        """Initialize the bounds for each parameter"""

        self.logger.info("Initializing bounds...\n")
        # We set the bounds for x0 and mu
        bounds = [
            self.conc_biom_bounds,
            self.flux_biom_bounds
        ]
        # If lag phase time should be estimated, add default bounds
        if self.t_lag:
            bounds.append(
                (0, 0.5 * self.time_vector.max())
            )
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
        self.logger.info(f"Bounds: {self.bounds}\n")

    def _build_weight_matrix(self):
        """
        Build the weight matrix from different input types

        :return: None
        """

        # First condition: the weights are in a 1D array
        if isinstance(self.sd, np.ndarray):
            # We first check that the weight vector is as long as the experimental matrix on the row axis
            if self.sd.size != self.experimental_matrix[0].size:
                raise ValueError("Weight vector not of right size")
            else:
                # We duplicate the vector column-wise to build a matrix of duplicated weight vectors
                self.sd = np.tile(self.sd, (len(self.experimental_matrix), 1))

        # Second condition: the weight is a scalar and must be broadcast to a matrix with same shape as the data
        elif isinstance(self.sd, int) or isinstance(self.sd, float):
            self.sd = np.full(self.experimental_matrix.shape, self.sd)
        else:
            raise RuntimeError("Unknown error")

    def _get_default_weights(self):
        """
        Build a default weight matrix. Default values:
            * Biomass: 0.2
            * Metabolites: 0.5
        :return: None
        """

        weights = [0.2]
        for name in range(len(self.name_vector) - 1):
            weights.append(0.5)
        self.sd = np.array(weights)
        self._build_weight_matrix()

    def optimize(self):
        """Run optimization and build the simulated matrix from the optimized parameters"""

        self.logger.info("\nRunning optimization...\n")
        self.optimize_results = PhysioFitter._run_optimization(self.params, self.simulate, self.experimental_matrix,
                                                               self.time_vector, self.deg_vector,
                                                               self.sd, self.bounds)
        self.parameter_stats = {
            "optimal": self.optimize_results.x
        }
        self.logger.info(f"Optimization results: \n{self.optimize_results}\n")
        for i, param in zip(self.ids, self.optimize_results.x):
            self.logger.info(f"\n{i} = {param}\n")
        self.simulated_matrix = self.simulate(self.optimize_results.x, self.experimental_matrix,
                                              self.time_vector, self.deg_vector)
        nan_sim_mat = np.copy(self.simulated_matrix)
        nan_sim_mat[np.isnan(self.experimental_matrix)] = np.nan
        self.simulated_data = DataFrame(data=nan_sim_mat, index=self.time_vector, columns=self.name_vector)
        self.simulated_data.index.name = "Time"
        self.logger.info(f"Final Simulated Data: \n{self.simulated_data}\n")

    @staticmethod
    def _simple_sim(params, exp_data_matrix, time_vector, deg):
        """Function to simulate the matrix using input parameters and the no lag analytical equation"""

        simulated_matrix = np.empty_like(exp_data_matrix)
        x_0 = params[0]
        mu = params[1]
        exp_mu_t = np.exp(mu * time_vector)
        simulated_matrix[:, 0] = x_0 * exp_mu_t
        for idx, i in enumerate(range(1, int(len(params) / 2))):
            q = params[i * 2]
            m_0 = params[i * 2 + 1]
            k = deg[idx]
            simulated_matrix[:, i] = q * (x_0 / (mu + k)) * (
                    np.exp(mu * time_vector) - np.exp(-k * time_vector)) + m_0 * np.exp(-k * time_vector)
        return simulated_matrix

    @staticmethod
    def _lag_sim(params, exp_data_matrix, time_vector, deg):
        """
        Function to simulate the matrix using input parameters and the no deg analytical equation (with lag)
        :param params:
        :param exp_data_matrix:
        :param time_vector: vector containing times
        :return:
        """

        simulated_matrix = np.empty_like(exp_data_matrix)
        x_0 = params[0]
        mu = params[1]
        t_lag = params[2]

        # We get indices in time vector where time < t_lag
        idx = np.nonzero(time_vector < t_lag)

        # Fill at those indices with x_0
        x_t_lag = np.full((len(idx) - 1,), x_0)

        # The rest of the biomass points are calculated as usual
        mult_by_time = x_0 * np.exp(mu * (time_vector[len(idx) - 1:] - t_lag))

        # Concatenate both vectors and transfer to X_t column of the simulated matrix
        simulated_matrix[:, 0] = np.concatenate((x_t_lag, mult_by_time), axis=None)
        exp_mu_t_lag = np.exp(mu * (time_vector - t_lag)) - 1

        for i in range(1, int(len(params) / 2)):
            q = params[i * 2 + 1]
            m_0 = params[i * 2 + 2]
            m_t_lag = np.full((len(idx) - 1,), m_0)
            mult_by_time = q * (x_0 / mu) * exp_mu_t_lag + m_0
            simulated_matrix[:, i] = np.concatenate((m_t_lag, mult_by_time), axis=None)
        return simulated_matrix

    @staticmethod
    def _total_sim(params, exp_data_matrix, time_vector, deg):
        """
        Function to simulate the matrix using input parameters and accounting for a lag phase and metabolite degradation

        :param params:
        :param exp_data_matrix:
        :param time_vector:
        :param deg:
        :return:
        """

        # Get end shape
        simulated_matrix = np.empty_like(exp_data_matrix)

        # Get initial params
        x_0 = params[0]
        mu = params[1]
        t_lag = params[2]

        # We get indices in time vector where time < t_lag
        idx = np.nonzero(time_vector < t_lag)

        # Fill at those indices with x_0
        x_t_lag = np.full((len(idx) - 1,), x_0)

        # The rest of the biomass points are calculated as usual
        mult_by_time = x_0 * np.exp(mu * (time_vector[len(idx) - 1:] - t_lag))
        simulated_matrix[:, 0] = np.concatenate((x_t_lag, mult_by_time), axis=None)

        for i in range(1, int(len(params) / 2)):
            q = params[i * 2 + 1]
            m_0 = params[i * 2 + 2]
            k = deg[idx]
            m_t_lag = np.full((len(idx) - 1,), m_0)
            mult_by_time = q * (x_0 / (mu + k)) * (
                    np.exp(mu * (time_vector - t_lag)) - np.exp(-k * (time_vector - t_lag))
            ) + (m_0 * np.exp(-k * time_vector))
            simulated_matrix[:, i] = np.concatenate((m_t_lag, mult_by_time), axis=None)
        return simulated_matrix

    @staticmethod
    def _calculate_cost(params, func, exp_data_matrix, time_vector, deg, weight_matrix):
        """Calculate the cost (residue) using the square of simulated-experimental over the SDs"""

        simulated_matrix = func(params, exp_data_matrix, time_vector, deg)
        cost_val = np.square((simulated_matrix - exp_data_matrix) / weight_matrix)
        residuum = np.nansum(cost_val)
        return residuum

    @staticmethod
    def _run_optimization(params, func, exp_data_matrix, time_vector, deg, weight_matrix, bounds):
        """
        Run the optimization on input parameters using the cost function and Scipy minimize (L-BFGS-B method
        that is deterministic and uses the gradient method for optimizing)
        """

        optimize_results = minimize(PhysioFitter._calculate_cost, x0=params, args=(
            func, exp_data_matrix, time_vector, deg, weight_matrix), method="L-BFGS-B", bounds=bounds,
                                    options={'maxcor': 30}
                                    )
        return optimize_results

    def monte_carlo_analysis(self):
        """
        Run a monte carlo analysis to calculate optimization standard deviations on parameters and simulated data points
        """

        if not self.optimize_results:
            raise RuntimeError("Running Monte Carlo simulation without having run the optimization is impossible "
                               "as best fit results are needed to generate the initial simulated matrix")

        self.logger.info(f"Running monte carlo analysis. Number of iterations: {self.iterations}\n")

        # Store the optimized results in variable that will be overridden on every pass
        opt_res = self.optimize_results
        opt_params_list = []
        matrices = []

        for _ in range(self.iterations):
            new_matrix = self._apply_noise()

            # We optimise the parameters using the noisy matrix as input
            opt_res = PhysioFitter._run_optimization(opt_res.x, self.simulate, new_matrix, self.time_vector,
                                                     self.deg_vector, self.sd, self.bounds)

            # Store the new simulated matrix in list for later use
            matrices.append(self.simulate(opt_res.x, new_matrix, self.time_vector, self.deg_vector))

            # Store the new optimised parameters in list for later use
            opt_params_list.append(opt_res.x)

        # Build a 3D array from all the simulated matrices to get standard deviation on each data point
        matrices = np.array(matrices)
        self.matrices_ci = {
            "lower_ci": np.percentile(matrices, 2.5, axis=0),
            "higher_ci": np.percentile(matrices, 97.5, axis=0)
        }

        # Compute the statistics on the list of parameters: means, sds, medians and confidence interval
        self._compute_parameter_stats(opt_params_list)
        self.logger.info(f"Optimized parameters statistics:")
        for key, value in self.parameter_stats.items():
            self.logger.info(f"{key}: {value}")

        # Apply nan mask to be coherent with the experimental matrix
        nan_lower_ci = np.copy(self.matrices_ci['lower_ci'])
        nan_higher_ci = np.copy(self.matrices_ci['higher_ci'])
        nan_lower_ci[np.isnan(self.experimental_matrix)] = np.nan
        nan_higher_ci[np.isnan(self.experimental_matrix)] = np.nan

        self.logger.info(f"Simulated matrix lower confidence interval:\n{nan_lower_ci}\n")
        self.logger.info(f"Simulated matrix higher confidence interval:\n{nan_higher_ci}\n")
        return

    def _compute_parameter_stats(self, opt_params_list):
        """
        Compute statistics on the optimized parameters from the monte carlo analysis.

        :param opt_params_list: list of optimized parameter arrays generated during the monte carlo analysis
        :return: parameter stats attribute containing means, sds, medians, low and high CI
        """

        opt_params_means = np.mean(np.array(opt_params_list), 0)
        opt_params_sds = np.std(np.array(opt_params_list), 0)
        opt_params_meds = np.median(np.array(opt_params_list), 0)
        conf_ints = np.column_stack((
            np.percentile(opt_params_list, 2.5, 0),
            np.percentile(opt_params_list, 97.5, 0)
        ))

        self.parameter_stats.update({
            "mean": opt_params_means,
            "sd": opt_params_sds,
            "median": opt_params_meds,
            "CI_2.5": conf_ints[:, 0],
            "CI_97.5": conf_ints[:, 1]
        })

        # self.parameter_stats_df = DataFrame()

    def khi2_test(self):

        number_measurements = np.count_nonzero(~np.isnan(self.experimental_matrix))
        number_params = len(self.params)
        dof = number_measurements - number_params
        cost = self._calculate_cost(self.optimize_results.x, self.simulate, self.experimental_matrix, self.time_vector,
                                    self.deg_vector, self.sd)
        p_val = chi2.cdf(cost, dof)

        khi2_res = {
            "khi2_value": cost,
            "number_of_measurements": number_measurements,
            "number_of_params": number_params,
            "Degrees_of_freedom": dof,
            "p_val": p_val
        }
        self.khi2_res = DataFrame.from_dict(khi2_res, "index", columns=["Values"])

        self.logger.info(f"khi2 test results:\n"
                         f"khi2 value: {cost}\n"
                         f"Number of measurements: {number_measurements}\n"
                         f"Number of parameters to fit: {number_params}\n"
                         f"Degrees of freedom: {dof}\n"
                         f"p-value = {p_val}\n")

        if p_val < 0.95:
            self.logger.info(f"At level of 95% confidence, the model fits the data good enough with respect to the "
                             f"provided measurement SD. Value: "
                             f"{p_val}")

        else:
            self.logger.info(f"At level of 95% confidence, the model does not fit the data good enough with respect to "
                             f"the provided measurement SD. "
                             f"Value: {p_val}\n")

    @staticmethod
    def _add_noise(vector, sd):
        """
        Add random noise to a given array using input standard deviations.

        :param vector: input array on which to apply noise
        :type vector: class: numpy.ndarray
        :param sd: standard deviation to apply to the input array
        :type sd: class: numpy.ndarray
        :return: noisy ndarray
        """

        output = np.random.default_rng().normal(loc=vector, scale=sd, size=vector.size)
        return output

    def _apply_noise(self):
        """
        Apply noise to the simulated matrix obtained using optimized parameters. SDs are obtained from the weight
        matrix
        """

        new_matrix = np.array([
            PhysioFitter._add_noise(self.simulated_matrix[idx, :], sd)
            for idx, sd in enumerate(self.sd)
        ])
        return new_matrix
