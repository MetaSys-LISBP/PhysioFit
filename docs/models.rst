Models
=========

User-made models
*****************

Overview
---------

Since PhysioFit 3.0.0, it is possible to create and implement models into the tool using your own python3 code. This
section of documentation will explain how to write your first model, how to test your model and implement it
into your local installation. A section will also describe how to submit your code for integration onto the
`Workflow4Metabolomics <https://workflow4metabolomics.usegalaxy.fr/>`_ platform for use in fluxomics workflows for example.

# TODO: In simulate part explain that you can use analytical functions or numerical differentiation (ODE) and adding extra functions in simulate

Creating your first model
--------------------------

Build the template
^^^^^^^^^^^^^^^^^^

To implement user-made models, PhysioFit leverages Python's object model to create classes that inherit from an Abstract
Base Class and that handles all the heavy-lifting for implementation. This lets the user follow a simple set of rules
to get their model working in the software.

The model you implement will actually be a class situated in a separate module. Start by opening a text file
using your IDE (Integrated Development Environment) of choice, and enter the following structure in the file::

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

This will be the base template from which you will build your model.

Populating the template
^^^^^^^^^^^^^^^^^^^^^^^

The first attribute to add in your model's *__init__* method is the model name. This name can be anything you wish, but it is strongly advised
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

We can now try testing that our model can be initialized properly. We shall use the block at the end of the file for
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

If you now run the file, you shall have a standard output in your console that resembles the following:

.. image:: _static/models/standard_out1.jpeg
   :scale: 100%

The next step is to prepare the parameters for simulations. There are two types of parameters that can be implemented
in a model: **parameters to estimate** and **fixed parameters**.

.. _parameters_to_estimate:

Parameters to estimate
""""""""""""""""""""""

The parameters to estimate are the parameters that must be optimized, and thus that will need bounds and initial values
to be initialized. We start by adding them to the get_params method, and giving them default starting values::

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

The next step is now to give the parameters to estimate default bounds for the optimization process. The bounds are a
class of objects that handle the logic and checks. They are derived from the python dict base class, and as such
implement the same methods such as update for example. Here is an example of how to implement them: ::

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

.. warning:: The keys in the bounds and in the parameters to estimate dictionary must correspond!

.. _fixed_parameters:

Fixed parameters
""""""""""""""""

The fixed parameters are parameters that are given as constants in the model equations. For example, in the case of
steady-state models that account for non enymatic degradation (see :ref:`default_steady-state_models`.), we need to give
the unstable metabolite a constant that will define it's rate of degradation ::

    self.fixed_parameters = {"Degradation": {
            metabolite: 2 for metabolite in self.metabolites
            }
        }

The different fixed parameters are given in a dictionary of dictionaries, where the first level is the name of the
parameter itself (here degradation) and the second level contains the mapping of metabolite-value pairs that will be
the default values initialized on launch (here we give a default value of 0 for every metabolite for example). Each
key of the first level will be used to initialize a panel to configure the values for the metabolites given in the
second level.

Adding a simulation function
""""""""""""""""""""""""""""

Once the different parameter sections have been written, the next step is to implement the simulation function that
will be called on each iteration of the optimization process (see :ref:`optimization_process` for more details).
To do this, first write out the function definition and insert the following parameters: ::

    @staticmethod
    def simulate(
            params_opti: list,
            data_matrix: np.ndarray,
            time_vector: np.ndarray,
            params_non_opti: dict
    ):
        pass

As shown, the function accepts 4 different arguments:
    * *params_opti*: a list containing the values for each parameter to estimate **in the order of apparition in the
      associated parameters_to_estimate dictionary** (see :ref:`parameters_to_estimate`)
    * *data_matrix*: the numpy array containing the experimental data (or data with the same shape)
    * *time_vector*: the numpy array containing the time points
    * *params_non_opti*: a dictionary containing the fixed parameters (see :ref:`fixed_parameters`)

Next you can start writing the body of the function. It is highly suggested to unpack the values from the
list of parameters to estimate into variables that possess the name of the associated parameter in the dictionary. To
get the right shape for the simulated matrix, one can use the *empty_like* function from the numpy library: ::

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

The math explaining the above simulation function can be found :ref:`here <default_steady-state_models>` (equations
5 and 6).

The above example showcases the use of analytical functions to simulate the flux dynamics. It is also possible to use
numerical differentiation functions (ODE). This may need the implementation of additional functions into the simulate
function. This can be done within the body of the simulate function: ::

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

As we can see, the solver needs a calculate derivative function to be passed in as an argument. This function is thus
created within the body of the simulate function, before it's use in the solver. More information on the mathematics
behind this implementation can be found :ref:`here <default_dynamic_models>`.


Adding the model to the GUI
---------------------------

One a model has been designed, it is time to test it out in the GUI. To integrate your model into the GUI, it is as
simple as copying the .py file and pasting it in your installation "models" folder. You can get the path towards the
models folder by opening a python kernel in your dedicated environment and initializing an IoHandler ::

    from physiofit.base.io import IoHandler
    io_handler = IoHandler()
    print(io_handler.get_local_model_folder())

.. note:: The model file name must follow the naming convention "model_[model number].py". If the last model in the list
          is the "model_5.py", the next one should be named "model_6.py".
