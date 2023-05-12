"""
Test the creation and use of the PhysioFitter
"""

import pytest
import physiofit as phyfit
from pandas import read_csv


def test_physiofitter(data, sds):
    """
    Test that the model and PhysioFitter can be safely instanciated from IoHandler
    """

    io = phyfit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", data)
    assert isinstance(model, phyfit.models.base_model.Model)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=sds,
        debug_mode=True
    )
    assert isinstance(fitter, phyfit.base.fitter.PhysioFitter)

def test_simulation(data, sds):
    pass
