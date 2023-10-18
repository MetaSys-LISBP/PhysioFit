Models
=========

User-made models
*****************

Overview
---------

Since PhysioFit 3.0.0, users can create and implement their own models to calculate fluxes and other growth parameters for any biological system. This
section explains how to write your first model, how to test the model and how to implement it
on your PhysioFit instance. We also detail how to submit your code for integration onto the
`Workflow4Metabolomics <https://workflow4metabolomics.usegalaxy.fr/>`_ platform for use in tailor-made workflows.

Creating your first model
--------------------------

Building the template
^^^^^^^^^^^^^^^^^^

To implement user-made models, PhysioFit leverages Python's object model to create classes that inherit from an Abstract
Base Class and that handles all the heavy-lifting for implementation. A simple set of rules enables
users to use their model in PhysioFit.

The model must be a class located in a dedicated module. Start by opening a text file
using your IDE (Integrated Development Environment), and enter the following structure in the file::

    from physiofit.models.base_model import Model

    class ChildModel(Model):

        def __init__(self, data):
            pass

        def get_params(self):
            pass

        @staticmethod
        def simulate():
            pass

    if __name__ == "__main__":
        pass

This is the base template to build your model. Methods *get_params* (to initialize and return model parameters) and *simulate* (to simulate metabolite dynamics for a given set of parameters) are mandatory. Additional methods are allowed if needed (e.g. to carry out intermediary steps for the simulation).

Adding equations and parameters
^^^^^^^^^^^^^^^^^^^^^^^

The first attribute to add in your model's *__init__* method is the model name. We strongly advise
to choose a name that helps the user understand what the model is destined to simulate. You must also add two other
attributes: the parameters to estimate & the fixed parameters. Finally, you must also call the *super().init(data)*
method to inherit the logic from the base class ::

    from physiofit.models.base_model import Model

    class ChildModel(Model):

        def __init__(self, data):
            super().__init__(data)
            self.model_name = "Tutorial model"
            self.parameters_to_estimate = None
            self.fixed_parameters = None

        def get_params(self):
            pass

        @staticmethod
        def simulate():
            pass

    if __name__ == "__main__":
        pass

.. note:: If your model does not contain fixed parameters, you must still initialize the attribute as None. This is
          considered good practice.

We can now try testing that the model can be initialized properly. Use the block at the end of the file for
testing purposes. Here is an example of how you can test the model::

    if __name__ == "__main__":

        test_data = pd.DataFrame(
            {
                "time": [0, 1, 2, 3],
                "experiment": ["A", "A", "A", "A"],
                "X": [0.5, 0.8, 1.2, 1.8],
                "Glucose": [12, 11.6, 11, 10.2]
            }
        )

        model = ChildModel(data=test_data)
        print(model)

If you now run the file, you should have a standard output in your console that looks like:

.. image:: _static/models/standard_out1.jpeg
   :scale: 100%

The next step is to prepare the parameters for simulations. PhysioFit allows two types of parameters: **parameters to estimate** and **fixed parameters**.

.. _parameters_to_estimate:

Parameters to estimate
""""""""""""""""""""""

The parameters to estimate are the parameters that must be estimated by PhysioFit, and thus that require bounds and initial values
to be initialized. The list of parameters and their default (initial) values must be returned by the *get_params* method::

    from physiofit.models.base_model import Model

    class ChildModel(Model):

        def __init__(self, data):
            super().__init__(data)
            self.model_name = "Tutorial model"
            self.parameters_to_estimate = None
            self.fixed_parameters = None

        def get_params(self):

            # Parameters are given in a dictionnary, where the key is
            # the parameter name and the value is a number that will
            # be the initial value for the optimization process

            self.parameters_to_estimate = {
                "BM_0": 1,
                "growth_rate": 1
            }

            # Do the same for all metabolite parameters to estimate
            # using a for loop:

            for metabolite in self.metabolites:
                self.parameters_to_estimate.update(
                    {
                        f"{metabolite}_flux" : 1,
                        f"{metabolite}_init_value" : 1
                    }
                )

        @staticmethod
        def simulate():
            pass

.. note:: For a given model, the number of metabolites may vary depending on the experiment, hence the metabolite-dependent parameters can be automatically defined in this function (as illustrated here using a for loop).

The next step is to define the default bounds used for the optimization process (these bounds can be changed in the GUI). The bounds are a
class of objects that handle the logic and checks. They are derived from the python dict base class, and as such
implement the same methods (e.g. update). Here is an example of how to implement the bounds: ::

    from physiofit.models.base_model import Model

    class ChildModel(Model):

        def __init__(self, data):
            super().__init__(data)
            self.model_name = "Tutorial model"
            self.parameters_to_estimate = None
            self.fixed_parameters = None

        def get_params(self):

            # Parameters are given in a dictionnary, where the key is
            # the parameter name and the value is a number that will
            # be the initial value for the optimization process

            self.parameters_to_estimate = {
                "BM_0": 1,
                "growth_rate": 1
            }

            # Instantiate the bounds object

            self.bounds = Bounds(
                {
                    "BM_0": (1e-3, 10),
                    "growth_rate": (1e-3, 3)
                }
            )

            # Do the same for all metabolite parameters to estimate
            # using a for loop:

            for metabolite in self.metabolites:
                self.parameters_to_estimate.update(
                    {
                        f"{metabolite}_flux" : 1,
                        f"{metabolite}_init_value" : 1
                    }
                )

                # Append the default bounds to the bounds attribute
                self.bounds.update(
                    {
                        f"{metabolite}_flux": (-50, 50),
                        f"{metabolite}_init_value": (1e-6, 50)
                    }
                )

        @staticmethod
        def simulate():
            pass

.. warning:: The keys in the bounds and in the parameters to estimate dictionary must be the same!

.. _fixed_parameters:

Fixed parameters
""""""""""""""""

The fixed parameters are parameters that are given as constants in the model equations and are not estimated by PhysioFit. For example, in the case of
steady-state models that account for non enymatic degradation (see :ref:`default_steady-state_models`.), we need to give
the unstable metabolite(s) a degradation constant (measured independently) ::

    self.fixed_parameters = {"Degradation": {
            metabolite: 2 for metabolite in self.metabolites
            }
        }

The different fixed parameters are given in a dictionary of dictionaries, where the first level is the name of the
parameter itself (here degradation) and the second level contains the mapping of metabolite-value pairs that will be
the default values initialized (here we give a default value of 2 for every metabolite for example). Each
key of the first level is used to initialize a widget in the GUI, thus allowing users to change the corresponding values for the metabolites given in the
second level.

Simulation function
"""""""""""""""""""

Once the *get_params* method has been implemented, the next step is to implement the simulation function that
will be called on each iteration of the optimization process to simulate the metabolite dynamics that correspond to a given set of parameters (see :ref:`optimization_process` for more details).
To do this, first write out the function definition: ::

    @staticmethod
    def simulate(
            params_opti: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            params_non_opti: dict
    ):
        pass

As shown above, this function takes four arguments:
    * *params_opti*: list containing the values of each parameter to estimate **in the order of apparition in the
      associated parameters_to_estimate dictionary** (see :ref:`parameters_to_estimate`)
    * *data_matrix*: numpy array containing the experimental data (or data with the same shape)
    * *time_vector*: numpy array containing the time points
    * *params_non_opti*: dictionary containing the fixed parameters (see :ref:`fixed_parameters`)

Now you can start writing the body of the function. For sake of clarity, we recommend unpacking parameters values from the 
list of parameters to estimate into internal variables. Th function *simulate* must return a matrix containing the simulation results, with the same shape as 
the matrix containing the experimental data. To initialize the simulated matrix, you can 
use the *empty_like* function from the numpy library: ::

    @staticmethod
    def simulate(
            params_opti: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            params_non_opti: dict
    ):
        # Get end shape
        simulated_matrix = np.empty_like(data_matrix)

        # Get initial params
        x_0 = params_opti[0]
        mu = params_opti[1]

        # Get X_0 values
        exp_mu_t = np.exp(mu * time_vector)
        simulated_matrix[:, 0] = x_0 * exp_mu_t
        fixed_params = [value for value in params_non_opti["Degradation"].values()]

        # Get parameter names and run the calculations column by column
        for i in range(1, int(len(params_opti) / 2)):
            q = params_opti[i * 2]
            m_0 = params_opti[i * 2 + 1]
            k = fixed_params[i - 1]
            exp_k_t = np.exp(-k * time_vector)
            simulated_matrix[:, i] = q * (x_0 / (mu + k)) \
                                     * (exp_mu_t - exp_k_t) \
                                     + m_0 * exp_k_t

        return simulated_matrix

The math corresponding to the simulation function provided above as example can be found :ref:`here <default_steady-state_models>` (equations
5 and 6).

This example showcases the use of analytical functions to simulate the flux dynamics. It is also possible to use
numerical derivation to solve a system of ordinary differential equations (ODEs), which can be usefull when algebric derivation is not straightforward. This require the implementation of additional functions into the simulate
function. The system of ODEs can be provided directly within the body of the simulate function: ::

    from scipy.integrate import solve_ivp

    @staticmethod
    def simulate(
            params_opti: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            params_non_opti: dict
    ):

        # Get parameters
        x_0 = params_opti[0]
        y_BM = params_opti[1]
        km = params_opti[2]
        qsmax = params_opti[3]
        s_0 = params_opti[4]
        y_P = params_opti[5]
        p_0 = params_opti[6]
        params = (y_BM, y_P, km, qsmax)

        # initialize variables at t=0
        state = [x_0, s_0, p_0]

        def calculate_derivative(t, state, y_BM, y_P, km, qsmax):

            # get substrate and biomass concentrations
            s_t = state[0]
            x_t = state[1]

            # calculate fluxes
            qs_t = qsmax * (s_t / (km + s_t))
            mu_t = y_BM * qs_t
            qp_t = y_P * qs_t

            # calculate derivatives
            dx = mu_t * x_t
            ds = -qs_t * x_t
            dp = qp_t * x_t

            return dx, ds, dp

        # simulate time-course concentrations
        sol = solve_ivp(
            fun=calculate_derivative,
            t_span=(np.min(time_vector), np.max(time_vector)),
            y0 = state,
            args=params,
            method="LSODA",
            t_eval = list(time_vector)
        )

        return sol.y.T

As we can see, the function *calculate_derivative* returns the derivatives of each metabolite concentration and is used by an ODEs solver that performs the simulations. This function is thus
created within the body of the simulate function, before being called by the solver. More information on the mathematics
behind this implementation can be found :ref:`here <default_dynamic_models>`.


Testing the model
---------------------------

One a model has been designed, it is time to test it! To integrate your model into the GUI, just copy the .py file 
in the folder 'models' of PhysioFit. You can get the path towards the
models folder by opening a python kernel in your dedicated environment and initializing an IoHandler ::

    from physiofit.base.io import IoHandler
    io_handler = IoHandler()
    print(io_handler.get_local_model_folder())

.. note:: The model file name must follow the naming convention "model_[model number].py". If the last model in the list
          is the "model_5.py", the next one should be named "model_6.py".

You can now launch PhysioFit's GUI, load a data file corresponding to the new model, select the model, and run flux calculation. In case of errors, 
have a look to the error message and correct the code.