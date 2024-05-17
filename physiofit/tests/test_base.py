"""
Test the creation and use of the PhysioFitter
"""
import logging

import numpy as np
import pandas as pd
import pytest
import physiofit

logging.getLogger("physiofit").setLevel(logging.ERROR)

def test_physiofitter(base_test_data):
    """
    Test that the model and PhysioFitter can be safely instantiated from
    IoHandler
    """

    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", base_test_data)
    assert isinstance(model, physiofit.models.base_model.Model)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=0.2,
        debug_mode=False
    )
    assert isinstance(fitter, physiofit.base.fitter.PhysioFitter)


def test_simulation(base_test_data):
    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", base_test_data)
    model.get_params()
    sim_mat = model.simulate(
        [param for param in model.parameters.values()],
        model.experimental_matrix,
        model.time_vector,
        model.args
    )
    assert isinstance(sim_mat, np.ndarray)


def test_wrong_entry_for_model_data(base_test_data):
    with pytest.raises(AttributeError):
        io = physiofit.base.io.IoHandler()
        model = io.select_model("Steady-state batch model", None)
    with pytest.raises(AttributeError):
        model = io.select_model(
            "Steady-state batch model",
            np.array(
                [[0, 1, 2], [0, 1, 2]]
            )
        )
    with pytest.raises(AttributeError):
        model = io.select_model(
            "Steady-state batch model",
            "Hello world this is an error"
        )


def test_optimization_process(base_test_data):
    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", base_test_data)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=0.2,
        debug_mode=False
    )
    fitter.optimize()
    assert isinstance(fitter.simulated_data, pd.DataFrame)
    for col in ["X", "Glucose", "Acetate"]:
        assert col in fitter.simulated_data.columns


def test_monte_carlo(base_test_data):
    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", base_test_data)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=0.2,
        debug_mode=False
    )
    fitter.optimize()
    fitter.monte_carlo_analysis()
    assert hasattr(fitter, "matrices_ci")
    assert isinstance(fitter.matrices_ci["lower_ci"], np.ndarray)
    assert isinstance(fitter.matrices_ci["higher_ci"], np.ndarray)
    assert np.subtract(fitter.matrices_ci["higher_ci"],
                       fitter.matrices_ci["lower_ci"]).all() > 0
    assert hasattr(fitter, "parameter_stats")
    assert isinstance(fitter.parameter_stats, dict)


def test_that_simulated_and_experimental_matrices_are_close(base_test_data):
    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", base_test_data)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=0.2,
        debug_mode=False
    )
    fitter.optimize()
    print(base_test_data)
    print(fitter.simulated_data)
    assert np.allclose(
        a=fitter.experimental_matrix,
        b=fitter.simulated_matrix,
        rtol=1e-3,
        equal_nan=True
    )
