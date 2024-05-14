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
            t_lag=(0, 0.5 * self.time_vector.max())
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

        self.args = {"Degradation": {
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

        # Get initial params
        x_0, mu, t_lag = parameters[:3]

        # We get indices in time vector where time < t_lag
        idx = np.nonzero(time_vector < t_lag)[0]

        # Fill at those indices with x_0
        x_t_lag = np.full(shape=(len(idx)), fill_value=x_0)

        # Get time vector from length of t_lag to end
        time_vector = time_vector[len(idx):]

        # The rest of the biomass points are calculated as usual
        mult_by_time = x_0 * np.exp(mu * (time_vector - t_lag))

        # Concatenate both vectors and transfer to X_t column of the
        # simulated matrix
        simulated_matrix[:, 0] = np.concatenate((x_t_lag, mult_by_time))

        # Get extra arguments
        arg_values = [value for value in args["Degradation"].values()]

        # optimize some calculations
        exp_mu_t_lag = np.exp(mu * (time_vector - t_lag))

        for i in range(1, len(parameters) // 2):
            q, m_0 = parameters[i * 2 + 1:i * 2 + 3]
            k = arg_values[i - 1]
            m_t_lag = np.full(len(idx), fill_value=m_0)
            mult_by_time = q * (x_0 / (mu + k)) * (exp_mu_t_lag - np.exp(
                        -k * time_vector)) + m_0 * np.exp(-k * time_vector)
            simulated_matrix[:, i] = np.concatenate((m_t_lag, mult_by_time))

        return np.clip(simulated_matrix, a_min=0, a_max=None)
