"""
Module containing the methods used by PhysioFit.
"""
from __future__ import annotations

import numpy as np

from physiofit.models.base_model import Model, Bounds


class ChildModel(Model):

    def __init__(self, data):

        super().__init__(data)
        self.name = "Steady-state batch model with lag phase and " \
                    "degradation of metabolites"
        self.vini = 1
        self.parameters = None
        self.args = None

    def get_params(self):

        self.parameters = {
            "X_0": self.vini,
            "growth_rate": self.vini,
            "t_lag": self.vini
        }
        self.bounds = Bounds(
            X_0=(1e-3, 10),
            growth_rate=(1e-3, 3),
            t_lag=(1e-6, 0.5 * self.time_vector.max())
        )
        for metabolite in self.metabolites:
            self.parameters.update(
                {
                    f"{metabolite}_q": self.vini,
                    f"{metabolite}_M0": self.vini
                }
            )
            self.bounds.update(
                {
                    f"{metabolite}_q": (-50, 50),
                    f"{metabolite}_M0": (1e-6, 50)
                }
            )

        self.args = {"Degradation constants": {
            met: 0 for met in self.metabolites
        }
        }

    @staticmethod
    def simulate(
            parameters: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            args: dict
    ):
        # Get end shape
        simulated_matrix = np.empty_like(data_matrix)
        # print("simulated_matrix", simulated_matrix)
        # print("simulated_matrix.shape", simulated_matrix.shape)

        # Get initial params
        x_0, mu, t_lag = parameters[:3]

        # We get indices in time vector where time < t_lag
        mask = (time_vector < t_lag)
        # print("mask", mask)

        # Get time vector from length of t_lag to end
        time_vector_lag = time_vector[mask]
        # print("len(time_vector_lag)", len(time_vector_lag))
        # print("time_vector_lag", time_vector_lag)
        time_vector_growth = time_vector[np.logical_not(mask)]
        # print("len(time_vector_growth)", len(time_vector_growth))
        # print("time_vector_growth", time_vector_growth)
        time_vector_diff = time_vector_growth - t_lag

        # optimize some calculations
        exp_mu_t_lag = np.exp(mu * time_vector_diff)

        # Fill at those indices with x_0
        x_t_lag = np.repeat(x_0, len(time_vector_lag))

        # The rest of the biomass points are calculated as usual
        mult_by_time = x_0 * exp_mu_t_lag

        # Concatenate both vectors and transfer to X_t column of the
        # simulated matrix
        simulated_matrix[:, 0] = np.concatenate((x_t_lag, mult_by_time))

        # Get extra arguments
        arg_values = [value for value in args[("Degradation "
                                               "constants")].values()]

        for i in range(1, len(parameters) // 2):
            q, m_0 = parameters[i * 2 + 1:i * 2 + 3]
            k = arg_values[i - 1]
            m_t_lag = m_0 * np.exp(-k * time_vector_lag)
            mult_by_time = q * (x_0 / (mu + k)) * (exp_mu_t_lag - np.exp(
                -k * time_vector_diff)) + m_0 * np.exp(-k *
                                                       time_vector_growth)
            simulated_matrix[:, i] = np.concatenate((m_t_lag, mult_by_time))

        return simulated_matrix
