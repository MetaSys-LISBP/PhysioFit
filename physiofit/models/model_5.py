"""
Module containing a dynamic growth and flux estimation model
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from physiofit.models.base_model import Model, Bounds


class ChildModel(Model):

    def __init__(self, data):

        super().__init__(data)
        self.name = "Dynamic Monod model (1 substrate, 1 product)"
        self.vini = 1
        self.parameters = None

    def get_params(self):

        self.parameters = {
            "X_0": self.vini,
            "y_BM": self.vini
        }

        self.bounds = Bounds(
            X_0=(1e-3, 1),
            y_BM=(1e-3, 1)
        )

        for metabolite in self.metabolites:
            if metabolite.startswith("S_"):
                self.parameters.update(
                    {
                        f"{metabolite}_km": self.vini,
                        f"{metabolite}_qsmax": self.vini,
                        f"{metabolite}_s_0": 100
                    }
                )
                self.bounds.update(
                    {
                        f"{metabolite}_km": (1e-6, 100),
                        f"{metabolite}_qsmax": (1e-6, 20),
                        f"{metabolite}_s_0": (1e-6, 100)
                    }
                )
                break

        for metabolite in self.metabolites:
            if metabolite.startswith("P_"):
                self.parameters.update(
                    {
                        f"{metabolite}_y_P": self.vini,
                        f"{metabolite}_p_0": 100
                    }

                )
                self.bounds.update(
                    {
                        f"{metabolite}_y_P": (1e-6, 2),
                        f"{metabolite}_p_0": (1e-6, 100)
                    }
                )
                break

        if len(self.parameters) != 7:
            raise ValueError(
                "This model expects 2 metabolites in the data file (1 "
                "substrate with name starting with 'S_' and 1 product with "
                "name starting with 'P_')."
            )

    @staticmethod
    def simulate(
            parameters: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            args: dict
    ):

        # Get parameters
        x_0, y_BM, km, qsmax, s_0, y_P, p_0 = parameters

        # initialize variables at t=0
        state = [x_0, s_0, p_0]

        def calculate_derivative(t, state, y_BM, y_P, km, qsmax):
            # get substrate and biomass concentrations
            x_t, s_t = state[0], state[1]

            # calculate fluxes
            qs_t = qsmax * (s_t / (km + s_t))
            mu_t = y_BM * qs_t
            qp_t = y_P * qs_t

            # calculate derivatives
            dx = mu_t * x_t
            ds = -qs_t * x_t
            dp = qp_t * x_t

            return dx, ds, dp

        # simulate time-course concentrations
        sol = solve_ivp(
            fun=calculate_derivative,
            t_span=(np.min(time_vector), np.max(time_vector)),
            y0=state,
            args=(y_BM, y_P, km, qsmax),
            method="LSODA",
            t_eval=list(time_vector)
        )

        return sol.y.T


if __name__ == "__main__":
    from physiofit.base.io import IoHandler

    io = IoHandler()
    data = io.read_data(
        r"C:\Users\legregam\PycharmProjects\PhysioFit\data\KEIO_test_data"
        r"\ode_test\KEIO_ROBOT6_1.tsv")
    data = data.sort_values("time")

    model = ChildModel(data)
    model.get_params()
    params = [param for param in model.parameters.values()]
    sol = ChildModel.simulate(
        params,
        model.data,
        model.time_vector,
        None
    )
    print(f"Times: {sol.t}")
    print(f"Values: {np.array(sol.y)}")
