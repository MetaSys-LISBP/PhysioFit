"""
PhysioFit software main module
"""

import logging

import numpy as np
from pandas import DataFrame
from scipy.optimize import minimize, differential_evolution
from scipy.stats import chi2

from physiofit.models.base_model import Model

logger = logging.getLogger(f"physiofit.{__name__}")



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

            residuum = sum((sim - meas) / sd)²

        * **optimization of the initial parameters** using `scipy.optimize.minimize ('Differential evolution', with polish with 'L-BFGS-B' method) <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html#optimize-minimize-lbfgsb>`_
        * **sensitivity analysis, khi2 tests and plotting**

    :param data: DataFrame containing data and passed by IOstream object
    :type data: class: pandas.DataFrame
    :param model: Model to initialize parameters and optimize
    :type model: physiofit.models.base_model.Model
    :param mc: Should Monte-Carlo sensitivity analysis be performed (default=True)
    :type mc: Boolean
    :param iterations: number of iterations for Monte-Carlo simulation (default=50)
    :type iterations: int
    :param sd: sd matrix used for residuum calculations. Can be:

                * a matrix with the same dimensions as the measurements matrix (but without the time column)
                * a named vector containing sds for all the metabolites provided in the input file
                * 0  in which case the matrix is automatically constructed from default values
                * a dictionary with the data column headers as keys and the associated value as a scalar or list

    :type sd: int, float, list, dict or ndarray
    """

    def __init__(
            self,
            data,
            model,
            mc=True,
            iterations=100,
            sd=None,
            debug_mode=False
    ):

        self.data = data
        self.model = model
        self.mc = mc
        self.iterations = iterations
        self.sd = sd
        self.debug_mode = debug_mode
        self.experimental_matrix = self.data.drop("time", axis=1).to_numpy()

        self.simulated_matrix = None
        self.simulated_data = None
        self.optimize_results = None
        self.simulate = None
        self.parameter_stats = None
        self.opt_params_sds = None
        self.matrices_ci = None
        self.opt_conf_ints = None
        self.khi2_res = None

    def verify_attrs(self):
        """Check that attributes are valid"""

        if type(self.iterations) is not int:
            raise TypeError(
                f"Number of monte carlo iterations must be an "
                f"integer, and not of type {type(self.iterations)}"
            )

        allowed_sds = [int, float, list, np.ndarray]
        if type(self.sd) not in allowed_sds:
            raise TypeError(
                f"sds is not in the right format ({type(self.sd)}. "
                f"Compatible formats are:\n{allowed_sds}"
            )

    def _sd_dict_to_matrix(self):
        """Convert sd dictionary to matrix/vector"""

        logger.debug("Initializing sd matrix from dict")
        # Perform checks
        for key in self.sd.keys():
            if key not in self.model.name_vector:
                raise KeyError(
                    f"The key {key} is not part of the data headers"
                )
        for name in self.model.name_vector:
            if name not in self.sd.keys():
                raise KeyError(
                    f"The key {name} is missing from the sds dict"
                )

        # Get lengths of each sd entry
        sd_lengths = [
            len(self.sd[key]) if type(self.sd[key]) not in [float, int] else 1
            for key in self.sd.keys()
        ]

        # Make sure that lengths are the same
        if not all(elem == sd_lengths[0] for elem in sd_lengths):
            raise ValueError("All sd vectors must have the same length")

        # Build matrix/vector
        if sd_lengths[0] == 1:
            self.sd = [self.sd[name] for name in self.model.name_vector]
        else:
            columns = (self.sd[name] for name in self.model.name_vector)
            matrix = np.column_stack(columns)
            self.sd = matrix

        logger.debug(f"SD object type: {type(self.sd)}")
        logger.debug(f"Array built from dict:\n{self.sd}")

    def initialize_sd_matrix(self):
        """
        Initialize the sd matrix from different types of inputs: single value,
        vector or matrix.

        :return: None
        """

        # This function can be optimized, if the input is a matrix we should
        # detect it directly
        logger.info("Initializing sd matrix...\n")

        # If sd is None, we generate the default matrix
        if self.sd is None or self.sd == {}:
            try:
                self.sd = {"X": 0.2}
                for col in self.data.columns[2:]:
                    self.sd.update({col: 0.5})
            except Exception:
                raise

        if isinstance(self.sd, dict):
            self._sd_dict_to_matrix()
        # When sd is a single value, we build a sd matrix containing the value
        # in all positions
        if isinstance(self.sd, int) or isinstance(self.sd, float):
            self._build_sd_matrix()
            logger.debug(f"SD matrix: {self.sd}\n")
            return
        if not isinstance(self.sd, np.ndarray) and not isinstance(self.sd,
                                                                  list):
            raise TypeError(
                f"Cannot coerce SD to array. Please check that a list or array "
                f"is given as input.\nCurrent input: \n{self.sd}"
            )
        else:
            self.sd = np.array(self.sd)
        if not np.issubdtype(self.sd.dtype, np.number):
            try:
                self.sd = self.sd.astype(float)
            except ValueError:
                raise ValueError(
                    f"The sd vector/matrix contains values that are not "
                    f"numeric. \nCurrent sd vector/matrix: \n{self.sd}"
                )
            except Exception as e:
                raise RuntimeError(f"Unknown error: {e}")
        else:
            # If the array is not the right shape, we assume it is a vector
            # that needs to be tiled into a matrix

            logger.debug(f"SD matrix: {self.sd}\n")
            logger.debug(f"Type of SD matrix: {type(self.sd)}")

            if self.sd.shape != self.experimental_matrix.shape:
                try:
                    self._build_sd_matrix()
                except ValueError:
                    raise
                except RuntimeError:
                    raise
            else:
                logger.debug(f"sd matrix: {self.sd}\n")
                return
        logger.info(f"sd Matrix:\n{self.sd}\n")

    def _build_sd_matrix(self):
        """
        Build the sd matrix from different input types

        :return: None
        """

        # First condition: the sds are in a 1D array
        if isinstance(self.sd, np.ndarray):
            # We first check that the sd vector is as long as the
            # experimental matrix on the row axis
            if self.sd.size != self.experimental_matrix[0].size:
                raise ValueError("sd vector not of right size")
            else:
                # We duplicate the vector column-wise to build a matrix of
                # duplicated sd vectors
                self.sd = np.tile(self.sd, (len(self.experimental_matrix), 1))

        # Second condition: the sd is a scalar and must be broadcast to a
        # matrix with same shape as the data
        elif isinstance(self.sd, int) or isinstance(self.sd, float):
            self.sd = np.full(self.experimental_matrix.shape, self.sd)
        else:
            raise RuntimeError("Unknown error")

    def _get_default_sds(self):
        """
        Build a default sd matrix. Default values:
            * Biomass: 0.2
            * Metabolites: 0.5
        :return: None
        """

        sds = [0.2]
        for name in range(len(self.model.name_vector) - 1):
            sds.append(0.5)
        self.sd = np.array(sds)
        self._build_sd_matrix()

    def optimize(self):
        """
        Run optimization and build the simulated matrix
        from the optimized parameters
        """

        logger.info("\nRunning optimization...\n")
        bounds = self.model.bounds()
        parameters = [param for param in self.model.parameters_to_estimate.values()]
        logger.debug(f"Simulate function = {self.model.simulate}")
        self.optimize_results = self._run_optimization(
            params=parameters,
            func=self.model.simulate,
            exp_data_matrix=self.experimental_matrix,
            time_vector=self.model.time_vector,
            non_opt_params=self.model.fixed_parameters,
            sd_matrix=self.sd,
            bounds=bounds,
            method="differential_evolution"
        )
        self.parameter_stats = {
            "optimal": self.optimize_results.x
        }
        logger.info(f"Optimization results: \n{self.optimize_results}\n")
        for i, param in zip(
                self.model.parameters_to_estimate, self.optimize_results.x
        ):
            logger.info(f"\n{i} = {param}\n")
        self.simulated_matrix = self.model.simulate(
            self.optimize_results.x,
            self.experimental_matrix,
            self.model.time_vector,
            self.model.fixed_parameters
        )
        logger.debug(f"simulated_matrix:\n{self.simulated_matrix}")
        nan_sim_mat = np.copy(self.simulated_matrix)
        nan_sim_mat[np.isnan(self.experimental_matrix)] = np.nan
        self.simulated_data = DataFrame(
            data=nan_sim_mat,
            index=self.model.time_vector,
            columns=self.model.name_vector
        )
        self.simulated_data.index.name = "Time"
        logger.info(f"Final Simulated Data: \n{self.simulated_data}\n")

    @staticmethod
    def _calculate_cost(
            params, func, exp_data_matrix, time_vector, non_opt_params,
            sd_matrix
    ):
        """
        Calculate the cost (residue) using the square of
        simulated-experimental over the SDs
        """

        simulated_matrix = func(params, exp_data_matrix, time_vector,
                                non_opt_params)
        cost_val = np.square((simulated_matrix - exp_data_matrix) / sd_matrix)
        residuum = np.nansum(cost_val)
        return residuum

    @staticmethod
    def _run_optimization(
            params: list,
            func: Model.simulate,
            exp_data_matrix: np.ndarray,
            time_vector: np.ndarray,
            non_opt_params: dict,
            sd_matrix: np.ndarray,
            bounds: tuple,
            method: str
    ):
        """
        Run the optimization on input parameters using the cost function and
        Scipy minimize (L-BFGS-B method that is deterministic and uses the
        gradient method for optimizing)
        """

        if method == "differential_evolution":
            optimize_results = differential_evolution(
                PhysioFitter._calculate_cost, bounds=bounds,
                args=(
                func, exp_data_matrix, time_vector, non_opt_params, sd_matrix),
                polish=True, x0=np.array(params)
            )
        elif method == "L-BFGS-B":
            optimize_results = minimize(
                PhysioFitter._calculate_cost, x0=np.array(params),
                args=(
                func, exp_data_matrix, time_vector, non_opt_params, sd_matrix),
                method="L-BFGS-B", bounds=bounds, options={'maxcor': 30}
            )
        else:
            raise ValueError(f"{method} is not implemented")

        return optimize_results

    def monte_carlo_analysis(self):
        """
        Run a monte carlo analysis to calculate optimization standard
        deviations on parameters and simulated data points
        """

        if not self.optimize_results:
            raise RuntimeError(
                "Running Monte Carlo simulation without having run the "
                "optimization is impossible as best fit results are needed to "
                "generate the initial simulated matrix"
            )

        logger.info(
            f"Running monte carlo analysis. Number of iterations: "
            f"{self.iterations}\n"
        )

        # Store the optimized results in variable that will be overridden on
        # every pass
        opt_res = self.optimize_results
        opt_params_list = []
        matrices = []

        for _ in range(self.iterations):
            new_matrix = self._apply_noise()

            # We optimise the parameters using the noisy matrix as input
            opt_res = PhysioFitter._run_optimization(
                opt_res.x, self.model.simulate, new_matrix, self.model.time_vector,
                self.model.fixed_parameters, self.sd, self.model.bounds(),
                "L-BFGS-B"
            )

            # Store the new simulated matrix in list for later use
            matrices.append(
                self.model.simulate(
                    opt_res.x, new_matrix, self.model.time_vector,
                    self.model.fixed_parameters
                )
            )

            # Store the new optimised parameters in list for later use
            opt_params_list.append(opt_res.x)

        # Build a 3D array from all the simulated matrices to get standard
        # deviation on each data point
        matrices = np.array(matrices)
        self.matrices_ci = {
            "lower_ci": np.percentile(matrices, 2.5, axis=0),
            "higher_ci": np.percentile(matrices, 97.5, axis=0)
        }

        # Compute the statistics on the list of parameters: means, sds,
        # medians and confidence interval
        self._compute_parameter_stats(opt_params_list)
        logger.info(f"Optimized parameters statistics:")
        for key, value in self.parameter_stats.items():
            logger.info(f"{key}: {value}")

        # Apply nan mask to be coherent with the experimental matrix
        nan_lower_ci = np.copy(self.matrices_ci['lower_ci'])
        nan_higher_ci = np.copy(self.matrices_ci['higher_ci'])
        nan_lower_ci[np.isnan(self.experimental_matrix)] = np.nan
        nan_higher_ci[np.isnan(self.experimental_matrix)] = np.nan

        logger.info(
            f"Simulated matrix lower confidence interval:\n{nan_lower_ci}\n"
        )
        logger.info(
            f"Simulated matrix higher confidence interval:\n{nan_higher_ci}\n"
        )
        return

    def _compute_parameter_stats(self, opt_params_list):
        """
        Compute statistics on the optimized parameters from the monte carlo
        analysis.

        :param opt_params_list: list of optimized parameter arrays generated
                                during the monte carlo analysis
        :return: parameter stats attribute containing means, sds, medians, low
                 and high CI
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

        number_measurements = np.count_nonzero(
            ~np.isnan(self.experimental_matrix)
        )
        number_params = len(self.model.parameters_to_estimate)
        dof = number_measurements - number_params
        cost = self._calculate_cost(
            self.optimize_results.x, self.model.simulate, self.experimental_matrix,
            self.model.time_vector, self.model.fixed_parameters, self.sd
        )
        p_val = chi2.cdf(cost, dof)

        khi2_res = {
            "khi2_value": cost,
            "number_of_measurements": number_measurements,
            "number_of_params": number_params,
            "Degrees_of_freedom": dof,
            "p_val": p_val
        }
        self.khi2_res = DataFrame.from_dict(
            khi2_res, "index", columns=["Values"]
        )

        logger.info(f"khi2 test results:\n"
                     f"khi2 value: {cost}\n"
                     f"Number of measurements: {number_measurements}\n"
                     f"Number of parameters to fit: {number_params}\n"
                     f"Degrees of freedom: {dof}\n"
                     f"p-value = {p_val}\n"
                    )

        if p_val < 0.95:
            logger.info(
                f"At level of 95% confidence, the model fits the data good "
                f"enough with respect to the provided measurement SD. "
                f"Value: {p_val}"
            )

        else:
            logger.info(
                f"At level of 95% confidence, the model does not fit the data "
                f"good enough with respect to the provided measurement SD. "
                f"Value: {p_val}\n"
            )

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

        output = np.random.default_rng().normal(
            loc=vector, scale=sd, size=vector.size
        )
        return output

    def _apply_noise(self):
        """
        Apply noise to the simulated matrix obtained using optimized
        parameters. SDs are obtained from the sd matrix
        """

        new_matrix = np.array([
            PhysioFitter._add_noise(self.simulated_matrix[idx, :], sd)
            for idx, sd in enumerate(self.sd)
        ])
        return new_matrix
