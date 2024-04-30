"""
This module contains tests for the models module.
"""
import numpy as np
from physiofit.base.io import IoHandler

def test_monod_model(monod_model_sds, monod_model_data):
    """
    Test that the Monod model using pyFOOMB simulated data & parameters given as input.
    """

    io = IoHandler()
    model = io.select_model("Dynamic Monod model (1 substrate, 1 product)", monod_model_data)
    model.get_params()
    fitter = io.initialize_fitter(
        data=model.data,
        model=model,
        sd=monod_model_sds,
        debug_mode=True
    )
    fitter.optimize()
    optimized_params = {
        name: param for name, param in zip(
            list(fitter.model.parameters_to_estimate.keys()),
            fitter.parameter_stats["optimal"]
        )
    }

    # X_0
    assert np.isclose(
        a=optimized_params["X_0"],
        b=0.01,
        atol=0.0005
    )
    # y_BM
    assert np.isclose(
        a=optimized_params["y_BM"],
        b=0.5,
        atol=0.025
    )
    # S_substrate_km
    assert np.isclose(
        a=optimized_params["S_substrate_km"],
        b=0.05,
        atol=0.0025
    )
    # S_substrate_qsmax
    assert np.isclose(
        a=optimized_params["S_substrate_qsmax"],
        b=0.8,
        atol=0.04
    )
    # S_substrate_s_0
    assert np.isclose(
        a=optimized_params["S_substrate_s_0"],
        b=20,
        atol=0.1
    )
    # P_product_y_P
    assert np.isclose(
        a=optimized_params["P_product_y_P"],
        b=0.2,
        atol=0.001
    )
