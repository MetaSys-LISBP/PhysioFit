"""
Module containing the methods used by PhysioFit.
"""

import numpy as np

from physiofit.models.base_model import Model


class Simple(Model):

    def __init__(self, data):

        super().__init__(data)
        self.method_name = "Simple Simulation"
        self.vini = 1

    def get_params(self):

        self.parameters_to_estimate = ["X_0", "mu"]
        self.bounds = {
            "X_0": (1e-3, 10),
            "mu": (1e-3, 3)
        }
        for metabolite in self.metabolites:
            self.parameters_to_estimate.append(f"{metabolite}_M0")
            self.parameters_to_estimate.append(f"{metabolite}_q")
            self.bounds.update(
                {
                    f"{metabolite}_M0": (1e-6, 50),
                    f"{metabolite}_q": (-50, 50)
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
        for idx, i in enumerate(range(1, int(len(params_opti) / 2))):
            q = params_opti[i * 2]
            m_0 = params_opti[i * 2 + 1]
            k = params_non_opti[idx]
            simulated_matrix[:, i] = q * (x_0 / (mu + k)) * (np.exp(
                mu * time_vector) - np.exp(-k * time_vector)) + m_0 * np.exp(
                -k * time_vector
            )

        return simulated_matrix
