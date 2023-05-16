"""
Test the creation and use of the PhysioFitter
"""
import numpy as np
import pandas as pd
import pytest
import physiofit


def test_physiofitter(data, sds):
    """
    Test that the model and PhysioFitter can be safely instanciated from IoHandler
    """

    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", data)
    assert isinstance(model, physiofit.models.base_model.Model)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=sds,
        debug_mode=True
    )
    assert isinstance(fitter, physiofit.base.fitter.PhysioFitter)

def test_simulation(data, sds):

    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", data)
    model.get_params()
    sim_mat = model.simulate(
        [param for param in model.parameters_to_estimate.values()],
        model.experimental_matrix,
        model.time_vector,
        model.fixed_parameters
    )
    assert isinstance(sim_mat, np.ndarray)

def test_wrong_entry_for_model_data(data, sds):

    with pytest.raises(AttributeError):
        io = physiofit.base.io.IoHandler()
        model = io.select_model("Steady-state batch model", None)
    with pytest.raises(AttributeError):
        model = io.select_model(
            "Steady-state batch model",
            np.array(
                [[0, 1, 2],[0, 1, 2]]
            )
        )
    with pytest.raises(AttributeError):
        model = io.select_model(
            "Steady-state batch model",
            "Hello world this is an error"
        )

def test_optimization_process(data, sds):

    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", data)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=sds,
        debug_mode=True
    )
    fitter.optimize()
    assert isinstance(fitter.simulated_data, pd.DataFrame)
    for col in ["X", "Glucose", "Acetate"]:
        assert col in fitter.simulated_data.columns

def test_monte_carlo(data, sds):
    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", data)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=sds,
        debug_mode=True
    )
    fitter.optimize()
    fitter.monte_carlo_analysis()
    assert hasattr(fitter, "matrices_ci")
    assert isinstance(fitter.matrices_ci["lower_ci"], np.ndarray)
    assert isinstance(fitter.matrices_ci["higher_ci"], np.ndarray)
    assert np.subtract(fitter.matrices_ci["higher_ci"], fitter.matrices_ci["lower_ci"]).all() > 0
    assert hasattr(fitter, "parameter_stats")
    assert isinstance(fitter.parameter_stats, dict)

def test_that_simulated_and_experimental_matrices_are_close(data, sds):
    io = physiofit.base.io.IoHandler()
    model = io.select_model("Steady-state batch model", data)
    model.get_params()
    fitter = io.initialize_fitter(
        model.data,
        model=model,
        sd=sds,
        debug_mode=True
    )
    fitter.optimize()
    assert np.allclose(
        a=fitter.experimental_matrix,
        b=fitter.simulated_data,
        rtol=1,
        equal_nan=True
    )

