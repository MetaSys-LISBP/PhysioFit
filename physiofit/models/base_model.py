from abc import ABC, abstractmethod

import pandas as pd
import numpy as np


# Initialization: nom des métabolites en entrée

# Get les parametres en fonction du modèle: retourne les noms des paramètres
# (optimisables et non),leurs valeurs initiales et les bounds pour ceux qui
# sont optimisables

# Besoin d'une fonction simulate propre au modèle qui prends les temps,
# les paramètres optimisables optimisables et non optimisables


class Model(ABC):

    def __init__(self, data: pd.DataFrame):
        self.data = data
        self.time_vector = self.data.time.to_numpy()
        self.name_vector = self.data.drop("time", axis=1).columns.to_list()
        self.experimental_matrix = self.data.drop("time", axis=1).to_numpy()
        self.metabolites = self.name_vector[1:]
        self.model_name = None
        self.parameters_to_estimate = None
        self.fixed_parameters = None
        self.bounds = None
        self.initial_values = None

    def __repr__(self):
        return f"Selected model: {self.model_name}\n" \
               f"Model data: \n{self.data}\n" \
               f"Experimental matrix: \n{self.experimental_matrix}\n" \
               f"Time vector: {self.time_vector}\n" \
               f"Name vector: {self.name_vector}\n" \
               f"Metabolites: {self.metabolites}\n" \
               f"Parameters to estimate: {self.parameters_to_estimate}\n" \
               f"Fixed parameters: {self.fixed_parameters}\n" \
               f"Bounds: {self.bounds}\n" \
               f"Initial values: {self.initial_values}\n" \
 \
    @abstractmethod
    def get_params(self):
        """

        :return params_to_estimate: List of parameters to estimate
        :return fixed_parameters: dict of constant parameters
        :return bounds: dict of upper and lower bounds
        :return default_init_values: dict containing default initial values for
                                    params
        """
        pass

    @staticmethod
    @abstractmethod
    def simulate(
            params_opti: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            params_non_opti: dict | list
    ):
        pass


class Bounds:

    __slots__ = "limits"
    def __init__(self, **kwargs):

        for key, value in kwargs.values():
            if not isinstance(value, tuple):
                raise TypeError(
                    "All bounds must be tuples"
                )
            for x in value:
                if not isinstance(x, int) or not isinstance(x, float):
                    raise TypeError(
                        "Bounds must be numbers"
                    )
            if value[0] >= value[1]:
                raise ValueError(
                    "Lower bound cannot be higher than upper bound. "
                )
            if not isinstance(key, str):
                raise TypeError(
                    "Name for bounds must be strings"
                )
        self.limits = kwargs

    def __repr__(self):
        for key, bound in self.limits.items():
            print(f"{key}: {bound}")

    def __getitem__(self, item):
        return self.limits[item]

    def __setitem__(self, key, value):
        self.limits[key] = value
