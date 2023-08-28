from __future__ import annotations

from abc import ABC, abstractmethod
from ast import literal_eval

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
        if "time" in self.data.columns and "experiments" in self.data.columns:
            self.name_vector = self.data.drop(["time", "experiments"], axis=1).columns.to_list()
        elif "time" in self.data.columns:
            self.name_vector = self.data.drop(["time"], axis=1).columns.to_list()
        else:
            raise ValueError(
                "Couldn't get column names from data. Is data loaded in properly and well formatted?"
            )
        self.experimental_matrix = self.data.drop("time", axis=1).to_numpy()
        self.metabolites = self.name_vector[1:]
        self.model_name = None
        self.parameters_to_estimate = None
        self.fixed_parameters = None
        self.bounds = None


    def __repr__(self):
        return f"Selected model: {self.model_name}\n" \
               f"Model data: \n{self.data}\n" \
               f"Experimental matrix: \n{self.experimental_matrix}\n" \
               f"Time vector: {self.time_vector}\n" \
               f"Name vector: {self.name_vector}\n" \
               f"Biomass & Metabolites: {self.metabolites}\n" \
               f"Parameters to estimate: {self.parameters_to_estimate}\n" \
               f"Fixed parameters: {self.fixed_parameters}\n" \
               f"Bounds: {self.bounds}\n"

    def __setattr__(self, key, value):
        if key == "data":
            self.__dict__["data"] = value
            self.time_vector = self.data.time.to_numpy()
            self.name_vector = self.data.drop("time", axis=1).columns.to_list()
            self.experimental_matrix = self.data.drop("time", axis=1).to_numpy()
            self.metabolites = self.name_vector[1:]
        else:
            self.__dict__[key] = value

    @ abstractmethod
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


class Bounds(dict):

    def __init__(self, mapping=None, **kwargs):

        if mapping is None:
            mapping = {}
        if kwargs:
            mapping.update(
                {key: value for key, value in kwargs.items()}
            )

        for key, value in mapping.items():
            key, value = self._check_bounds(key, value)
            mapping[key] = value

        super().__init__(mapping)

    def __call__(self):

        tuples = tuple(value for value in self.values())
        return tuples

    def __setitem__(self, key, value):

        key, value = self._check_bounds(key, value)
        super().__setitem__(key, value)

    @staticmethod
    def _check_bounds(key, value):

        if not isinstance(value, tuple):
            if isinstance(value, str):
                try:
                    value = literal_eval(value)
                except Exception:
                    raise TypeError(
                        f"Could not coerce {value} into string"
                    )
        for x in value:
            if not isinstance(x, int) and not isinstance(x, float):
                raise TypeError(
                    "Individual bound values must be numbers"
                )
        if value[0] >= value[1]:
            raise ValueError(
                "Lower bound cannot be higher than upper bound. "
            )
        if len(value) != 2:
            raise ValueError(
                "Bounds can only have two values. Number of values detected: "
                f"{len(value)}"
            )
        if not isinstance(key, str):
            raise TypeError(
                "Name for bounds must be strings"
            )
        return key, value


class StandardDevs(dict):

    def __init__(self, mapping=None, **kwargs):

        self._vector = None

        if mapping is None:
            mapping = {}
        if kwargs:
            mapping.update(
                {key: value for key, value in kwargs.items()}
            )

        for key, value in mapping.items():
            self._check_sd(key, value)

        super().__init__(mapping)

    @staticmethod
    def _check_sd(key, value):

        if not isinstance(key, str):
            raise TypeError(
                f"SD name field can only contain strings. Detected type {type(key)} for {key}"
            )
        if not isinstance(value, int) and not isinstance(value, float):
            try:
                value = float(value)
            except Exception:
                raise TypeError(
                    f"SD value must be a number. Detected type: {type(value)} for {key}"
                )
        if value <= 0:
            raise ValueError(
                f"SD value must be superior to 0. Detected value: {value} for {key}"
            )

    def __setitem__(self, key, value):

        self._check_sd(key, value)
        super().__setitem__(key, value)

    @property
    def vector(self):

        if self._vector is not None:
            return self._vector
        self._vector = np.array([value for value in self.values])
        return self._vector

class ModelError(Exception):
    pass


if __name__ == "__main__":
    bounds = Bounds(
        X_0=(1e-3, 10),
        mu=(1e-3, 3)
    )
    print(bounds())
