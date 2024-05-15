"""
Configuration variables for the tests

Shared pre-defined parameters
The return value of fixture function will be available as a predefined
parameter for all test functions. The test's parameter name must be the same as
the fixture function's name.
"""
import numpy as np
from numpy import array
import pandas as pd
import pytest

from physiofit.models.base_model import StandardDevs
from physiofit.base.io import IoHandler


@pytest.fixture
def placeholder_data():
    return pd.DataFrame(
        {
            "time": np.arange(6, step=0.2),
            "X": np.arange(6, step=0.2),
            "Glucose": np.arange(6, step=0.2),
            "Acetate": np.arange(6, step=0.2),
            "Glutamate": np.arange(6, step=0.2),
        }
    )


@pytest.fixture
def base_test_data():
    """Test data to use in tests
    (taken from Bergès et al., 2021 --> KEIO_ROBOT1_1)"""

    return pd.DataFrame.from_dict(
        {
            "time": [0, 1.18, 2.27, 3.13, 3.77, 4.42,
                     4.82, 0.67, 1.72, 2.8, 3.63, 4.27, 4.88],
            "X": [0.03, 0.05, 0.08, 0.13, 0.18, 0.24, 0.34,
                  np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            "Glucose": [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                        14, 15.30, 13.68, 12.81, 12.15, 10.93],
            "Acetate": [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                        0.01, 0.33, 0.72, 1.17, 1.63, 2.14]
        }
    )


@pytest.fixture
def parameters():
    return {
        "t_lag": 1.4,
        "X_0": 0.02,
        "growth_rate": 0.8,
        "Glucose_q": -8,
        "Glucose_M0": 20,
        "Acetate_q": 3,
        "Acetate_M0": 0.01,
        "Glutamate_q": 2,
        "Glutamate_M0": 0.01
    }


@pytest.fixture
def model_2_data(placeholder_data, parameters):
    """ Test data to use in tests for the model_2: Steady-state batch model
    with lag phase. Data is simulated using synthetic parameters"""

    io = IoHandler()
    model = io.select_model(
        name="Steady-state batch model with lag phase",
        data=placeholder_data
    )
    model.get_params()
    model.parameters.update(parameters)
    sim_data = model.simulate(
        list(model.parameters.values()),
        model.data.drop("time", axis=1),
        model.time_vector,
        model.args
    )
    df = pd.DataFrame(
        data=sim_data,
        index=model.time_vector,
        columns=model.name_vector
    )
    df.index.name = "time"
    return df.reset_index()


@pytest.fixture
def model_4_data(placeholder_data, parameters):
    """ Test data to use in tests for the model_2: Steady-state batch model
    with lag phase. Data is simulated using synthetic parameters"""

    io = IoHandler()
    model = io.select_model(
        name="Steady-state batch model",
        data=placeholder_data
    )
    model.get_params()
    del parameters["t_lag"]
    model.parameters.update(parameters)
    sim_data = model.simulate(
        list(model.parameters.values()),
        model.data.drop("time", axis=1),
        model.time_vector,
        model.args
    )
    df = pd.DataFrame(
        data=sim_data,
        index=model.time_vector,
        columns=model.name_vector
    )
    df.index.name = "time"
    return df.reset_index()


@pytest.fixture
def monod_model_sds():
    return StandardDevs(
        X=0.2,
        S_substrate=0.2,
        P_product=0.2
    )


@pytest.fixture
def sds():
    return StandardDevs(
        X=0.2,
        Glucose=0.2,
        Acetate=0.2,
        Glutamate=0.2
    )
