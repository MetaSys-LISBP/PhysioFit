
"""
Configuration variables for the tests

Shared pre-defined parameters
The return value of fixture function will be available as a predefined
parameter for all test functions. The test's parameter name must be the same as
the fixture function's name.
"""

import os

import pandas as pd
import pytest
from numpy import nan
from physiofit.base.io import IoHandler
from physiofit.models.base_model import StandardDevs

@pytest.fixture
def data():
    """Test data to use in tests (taken from Bergès et al., 2021 --> KEIO_ROBOT1_1)"""

    return pd.DataFrame.from_dict(
        {
            "time": [0, 1.18, 2.27, 3.13, 3.77, 4.42, 4.82, 0.67, 1.72, 2.8, 3.63, 4.27, 4.88],
            "X": [0.03, 0.05, 0.08, 0.13, 0.18, 0.24, 0.34, nan, nan, nan, nan, nan, nan],
            "Glucose": [nan, nan, nan, nan, nan, nan, nan, 14, 15.30, 13.68, 12.81, 12.15, 10.93],
            "Acetate": [nan, nan, nan, nan, nan, nan, nan, 0.01, 0.33, 0.72, 1.17, 1.63, 2.14]
        }
    )

# time	X	Glc	Ace
# 0	0.031752	NA	NA
# 1.18888888888889	0.05292	NA	NA
# 2.27694444444444	0.084672	NA	NA
# 3.12833333333333	0.134568	NA	NA
# 3.77138888888889	0.175392	NA	NA
# 4.41555555555556	0.244944	NA	NA
# 4.82277777777778	0.337176	NA	NA
# 0.0666666666666667	NA	14.0042771787617	0.0130592124021829
# 1.71666666666667	NA	15.3032563161238	0.325073443768164
# 2.8	NA	13.6751055268759	0.722570137805818
# 3.63333333333333	NA	12.8110643342243	1.16940283499877
# 4.26666666666667	NA	12.1520761714542	1.62777664140305
# 4.88333333333333	NA	10.9297935900362	2.13639255382759

@pytest.fixture
def sds():
    return StandardDevs(
        X = 0.2,
        Glucose = 0.5,
        Acetate = 0.5
    )

