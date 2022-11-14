"""
Module containing the methods used by PhysioFit.
"""

import numpy as np

from physiofit.models.base_model import Model
from ..base.io import Bounds


class ChildModel(Model):

    def __init__(self, data):

        super().__init__(data)
        self.model_name = "General model"
        self.vini = 1
        self.parameters_to_estimate = None
        self.initial_values = None

    def get_params(self):

        self.parameters_to_estimate = ["X_0", "mu"]
        self.bounds = {
            "X_0": Bounds(1e-3, 10),
            "mu": Bounds(1e-3, 3)
        }
        for metabolite in self.metabolites:
            self.parameters_to_estimate.append(f"{metabolite}_q")
            self.parameters_to_estimate.append(f"{metabolite}_M0")
            self.bounds.update(
                {
                    f"{metabolite}_q": Bounds(-50, 50),
                    f"{metabolite}_M0": Bounds(1e-6, 50)
                }
            )
        self.initial_values = {
            i: self.vini for i in self.parameters_to_estimate
        }

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
        exp_mu_t = np.exp(mu * time_vector)
        simulated_matrix[:, 0] = x_0 * exp_mu_t
        for i in range(1, int(len(params_opti) / 2)):
            q = params_opti[i * 2]
            m_0 = params_opti[i * 2 + 1]
            k = params_non_opti[i - 1]
            exp_k_t = np.exp(-k * time_vector)
            simulated_matrix[:, i] = q * (x_0 / (mu + k)) \
                                     * (exp_mu_t - exp_k_t) \
                                     + m_0 * exp_k_t

        return simulated_matrix
