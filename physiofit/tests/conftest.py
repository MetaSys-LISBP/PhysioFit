
"""
Configuration variables for the tests

Shared pre-defined parameters
The return value of fixture function will be available as a predefined
parameter for all test functions. The test's parameter name must be the same as
the fixture function's name.
"""

import pytest
from physiofit.base.io import IoHandler
from physiofit.models.base_model import StandardDevs

@pytest.fixture
def data():
    """Test data to use in tests"""

    return IoHandler.read_data(
        r"C:\Users\legregam\PycharmProjects\PhysioFit\physiofit\data\data_example.tsv"
    )

@pytest.fixture
def sds():
    return StandardDevs(
        X = 0.2,
        Glc = 0.5,
        Ace = 0.5
    )

