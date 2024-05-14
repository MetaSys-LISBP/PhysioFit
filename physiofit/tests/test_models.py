"""
This module contains tests for the models module.
"""
import numpy as np
import pandas as pd

from physiofit.base.io import IoHandler
import physiofit.models.base_model


def test_steady_state_batch_model_with_lag_phase_and_degradation_of_metabolites(
        model_1_data: pd.DataFrame,
        sds: physiofit.models.base_model.StandardDevs
):
    io = IoHandler()
    model = io.select_model(
        name="Steady-state batch model with lag phase and "
        "degradation of metabolites",
        data=model_1_data
    )
    model.get_params()
    fitter = io.initialize_fitter(
        data=model.data,
        model=model,
        sd=sds,
        debug_mode=True
    )
    fitter.optimize()
    optimized_params = {
        name: param for name, param in zip(
            list(fitter.model.parameters.keys()),
            fitter.parameter_stats["optimal"]
        )
    }

    # t_lag
    assert np.isclose(
        a=optimized_params["t_lag"],
        b=1.4,
        rtol=0.2
    )
    # X_0
    assert np.isclose(
        a=optimized_params["X_0"],
        b=0.02,
        rtol=0.02
    )
    # growth_rate
    assert np.isclose(
        a=optimized_params["growth_rate"],
        b=0.8,
        rtol=0.02
    )
    # Glucose_q
    assert np.isclose(
        a=optimized_params["Glucose_q"],
        b=-8,
        rtol=0.02
    )
    # Glucose_M0
    assert np.isclose(
        a=optimized_params["Glucose_M0"],
        b=20,
        rtol=0.02
    )
    # Acetate_q
    assert np.isclose(
        a=optimized_params["Acetate_q"],
        b=3,
        rtol=0.02
    )
    # Acetate_M0
    assert np.isclose(
        a=optimized_params["Acetate_M0"],
        b=0.01,
        rtol=0.02
    )
    # Glutamate_q
    assert np.isclose(
        a=optimized_params["Glutamate_q"],
        b=0.05,
        rtol=0.02
    )
    # Glutamate_M0
    assert np.isclose(
        a=optimized_params["Glutamate_M0"],
        b=0.01,
        rtol=0.02
    )


def test_steady_state_batch_model_with_lag_phase(
        model_2_data: pd.DataFrame,
        sds: physiofit.models.base_model.StandardDevs
):
    io = IoHandler()
    model = io.select_model(
        name="Steady-state batch model with lag phase",
        data=model_2_data
    )
    model.get_params()
    fitter = io.initialize_fitter(
        data=model.data,
        model=model,
        sd=sds,
        debug_mode=True
    )
    fitter.optimize()
    optimized_params = {
        name: param for name, param in zip(
            list(fitter.model.parameters.keys()),
            fitter.parameter_stats["optimal"]
        )
    }

    # t_lag
    assert np.isclose(
        a=optimized_params["t_lag"],
        b=1.4,
        rtol=0.01
    )
    # X_0
    assert np.isclose(
        a=optimized_params["X_0"],
        b=0.02,
        rtol=0.01
    )
    # growth_rate
    assert np.isclose(
        a=optimized_params["growth_rate"],
        b=0.8,
        rtol=0.01
    )
    # Glucose_q
    assert np.isclose(
        a=optimized_params["Glucose_q"],
        b=-8,
        rtol=0.01
    )
    # Glucose_M0
    assert np.isclose(
        a=optimized_params["Glucose_M0"],
        b=20,
        rtol=0.01
    )
    # Acetate_q
    assert np.isclose(
        a=optimized_params["Acetate_q"],
        b=3,
        rtol=0.01
    )
    # Acetate_M0
    assert np.isclose(
        a=optimized_params["Acetate_M0"],
        b=0.01,
        rtol=0.01
    )
    # Glutamate_q
    assert np.isclose(
        a=optimized_params["Glutamate_q"],
        b=2,
        rtol=0.01
    )
    # Glutamate_M0
    assert np.isclose(
        a=optimized_params["Glutamate_M0"],
        b=0.01,
        rtol=0.01
    )

def test_steady_state_batch_model(
        model_4_data: pd.DataFrame,
        sds: physiofit.models.base_model.StandardDevs
):
    io = IoHandler()
    model = io.select_model(
        "Steady-state batch model",
        model_4_data
    )
    model.get_params()
    fitter = io.initialize_fitter(
        data=model.data,
        model=model,
        sd=sds,
        debug_mode=True
    )
    fitter.optimize()
    optimized_params = {
        name: param for name, param in zip(
            list(fitter.model.parameters.keys()),
            fitter.parameter_stats["optimal"]
        )
    }
    # X_0
    assert np.isclose(
        a=optimized_params["X_0"],
        b=0.02,
        rtol=0.01
    )
    # growth_rate
    assert np.isclose(
        a=optimized_params["growth_rate"],
        b=0.8,
        rtol=0.01
    )
    # Glucose_q
    assert np.isclose(
        a=optimized_params["Glucose_q"],
        b=-8,
        rtol=0.01
    )
    # Glucose_M0
    assert np.isclose(
        a=optimized_params["Glucose_M0"],
        b=20,
        rtol=0.01
    )
    # Acetate_q
    assert np.isclose(
        a=optimized_params["Acetate_q"],
        b=3,
        rtol=0.01
    )
    # Acetate_M0
    assert np.isclose(
        a=optimized_params["Acetate_M0"],
        b=0.01,
        rtol=0.01
    )
    # Glutamate_q
    assert np.isclose(
        a=optimized_params["Glutamate_q"],
        b=2,
        rtol=0.01
    )
    # Glutamate_M0
    assert np.isclose(
        a=optimized_params["Glutamate_M0"],
        b=0.01,
        rtol=0.01
    )


def test_monod_model(monod_model_sds, monod_model_data):
    """
    Test that the Monod model using pyFOOMB simulated data & parameters given
    as input.
    """

    io = IoHandler()
    model = io.select_model(
        "Dynamic Monod model (1 substrate, 1 product)",
        monod_model_data
    )
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
            list(fitter.model.parameters.keys()),
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
        b=0.1,
        atol=0.001
    )
