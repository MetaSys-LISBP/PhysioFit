"""
Module containing the methods used by PhysioFit.
"""

from __future__ import annotations
import numpy as np

from physiofit.models.base_model import Model, Bounds


class ChildModel(Model):

    def __init__(self, data):

        super().__init__(data)
        self.name = "Steady-state batch model"
        self.vini = 1
        self.parameters = None
        self.args = None

    def get_params(self):

        self.parameters = {
            "X_0": self.vini,
            "growth_rate": self.vini
        }
        self.bounds = Bounds(
            {
                "X_0": (1e-3, 10),
                "growth_rate": (1e-3, 3)
            }
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

    @staticmethod
    def simulate(
            parameters: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            args: dict | list
    ):
        # Get end shape
        simulated_matrix = np.empty_like(data_matrix)

        # Get parameters
        x_0, mu = parameters[:2]

        # Build matrix
        exp_mu_t = np.exp(mu * time_vector)
        simulated_matrix[:, 0] = x_0 * exp_mu_t
        for i in range(1, len(parameters) // 2):
            q, m_0 = parameters[i * 2:i * 2 + 2]
            simulated_matrix[:, i] = q * (x_0 / mu) * (exp_mu_t - 1) + m_0

        return np.clip(simulated_matrix, a_min=0, a_max=None)
