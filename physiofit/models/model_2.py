"""
Module containing the methods used by PhysioFit.
"""

from __future__ import annotations
import numpy as np

from physiofit.models.base_model import Model, Bounds


class ChildModel(Model):

    def __init__(self, data):

        super().__init__(data)
        self.model_name = "Steady-state batch model with lag phase"
        self.vini = 1
        self.parameters_to_estimate = None
        self.fixed_parameters = None

    def get_params(self):

        self.parameters_to_estimate = {
            "X_0" : self.vini,
            "growth_rate" : self.vini,
            "t_lag" : self.vini
        }
        self.bounds = Bounds({
            "X_0": (1e-3, 10),
            "growth_rate": (1e-3, 3),
            "t_lag": (0, 0.5 * self.time_vector.max()),
        })
        for metabolite in self.metabolites:
            self.parameters_to_estimate.update(
                {
                    f"{metabolite}_q" : 1,
                    f"{metabolite}_M0" : 1
                }
            )
            self.bounds.update(
                {
                    f"{metabolite}_q": (-50, 50),
                    f"{metabolite}_M0": (1e-6, 50)
                }
            )

    @staticmethod
    def simulate(
            params_opti: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            params_non_opti: dict | list
    ):
        # Get end shape
        simulated_matrix = np.empty_like(data_matrix)
        x_0 = params_opti[0]
        mu = params_opti[1]
        t_lag = params_opti[2]

        # We get indices in time vector where time < t_lag
        idx = np.nonzero(time_vector < t_lag)

        # Fill at those indices with x_0
        x_t_lag = np.full((len(idx) - 1,), x_0)

        # The rest of the biomass points are calculated as usual
        mult_by_time = x_0 * np.exp(mu * (time_vector[len(idx) - 1:] - t_lag))

        # Concatenate both vectors and transfer to X_t column of the
        # simulated matrix
        simulated_matrix[:, 0] = np.concatenate(
            (x_t_lag, mult_by_time),
            axis=None
        )
        exp_mu_t_lag = np.exp(mu * (time_vector - t_lag)) - 1

        for i in range(1, int(len(params_opti) / 2)):
            q = params_opti[i * 2 + 1]
            m_0 = params_opti[i * 2 + 2]
            m_t_lag = np.full((len(idx) - 1,), m_0)
            mult_by_time = q * (x_0 / mu) * exp_mu_t_lag + m_0
            simulated_matrix[:, i] = np.concatenate(
                (m_t_lag, mult_by_time),
                axis=None
            )

        return simulated_matrix
