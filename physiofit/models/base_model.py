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
        self.params_to_estimate = None
        self.fixed_parameters = None
        self.bounds = None
        self.default_init_values = None

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

