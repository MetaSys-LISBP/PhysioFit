"""
Module containing the methods used by PhysioFit.
"""

from __future__ import annotations
import numpy as np

from physiofit.models.base_model import Model, Bounds


class ChildModel(Model):

    def __init__(self, data):

        super().__init__(data)
        self.name = ("Steady-state batch model with degradation of "
                     "metabolites")
        self.vini = 1
        self.parameters = None
        self.args = None

    def get_params(self):

        self.parameters = {
            "X_0": self.vini,
            "growth_rate": self.vini
        }
        self.args = {"Degradation constant": {
            met: 0 for met in self.metabolites
            }
        }
        self.bounds = Bounds({
            "X_0": (1e-3, 10),
            "growth_rate": (1e-3, 3),
        })
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
            args: dict | list
    ):
        # Get end shape
        simulated_matrix = np.empty_like(data_matrix)

        # Get initial params
        x_0 = parameters[0]
        mu = parameters[1]

        # Get X_0 values
        exp_mu_t = np.exp(mu * time_vector)
        simulated_matrix[:, 0] = x_0 * exp_mu_t
        fixed_params = [value for value in args[("Degradation "
                                                 "constants")].values()]

        for i in range(1, int(len(parameters) / 2)):
            q = parameters[i * 2]
            m_0 = parameters[i * 2 + 1]
            k = fixed_params[i - 1]
            exp_k_t = np.exp(-k * time_vector)
            simulated_matrix[:, i] = q * (x_0 / (mu + k)) \
                                     * (exp_mu_t - exp_k_t) \
                                     + m_0 * exp_k_t

        return simulated_matrix
