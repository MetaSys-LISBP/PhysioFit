Models
=========


Models shipped with PhysioFit
******************************

.. _default_steady-state_models:

Steady-state models
----------------------------

As detailed in Peiro et al. (2019), PhysioFit includes a steady-state model that accounts for i) non enzymatic degradation of some metabolites and
ii) growth lag. This model is described by the following system of ordinary differential equations:

.. math:: 
  \begin{align}
    \dfrac{dX}{dt} = \begin{cases}
    0 & \text{if } t < t_{lag} \\
    \mu \cdot X & \text{ow.}
  \end{cases}    \label{eq1}    \tag{eq. 1} 
  \end{align}

where :math:`\mu` is growth rate, :math:`X` is the biomass concentration, and :math:`t_{lag}` is the lag time. 

.. math:: 
  \begin{align}
    \dfrac{dM_{i}}{dt} = \begin{cases}
    -k_{i}\cdot M_{i} & \text{if } t < t_{lag} \\
    -k_{i}\cdot M_{i}+X\cdot qM_{i} & \text{ow.}
  \end{cases}    \label{eq2}    \tag{eq. 2} 
  \end{align}

where :math:`k_{i}` is the first-order degradation constant of the 
metabolite :math:`M_{i}` and :math:`qM_{i}` is its exchange (uptake or production) flux. 

The flux :math:`qM_{i}` is positive (negative) when :math:`M_{i}` is produced (consumed). The sign of :math:`qM_{i}` can thus be used to
automatically identify products and substrates.

Integrating equations 1-2 provides the following analytical functions:

.. math:: 
  \begin{align}
    X(t) = \begin{cases}
    X_{0} & \text{if } t < t_{lag} \\
    X_{0} \cdot e^{\mu \cdot (t - t_{lag})} & \text{ow.}
  \end{cases}    \label{eq3}    \tag{eq. 3} 
  \end{align}

.. math:: 
  \begin{align}
    M_{i}(t) = \begin{cases}
    M_{i}^{0} \cdot e^{-k_{i} \cdot t} & \text{if } t < t_{lag} \\
    qM_{i} \cdot \dfrac{X_{0}}{\mu + k_{i}} \cdot (e^{\mu \cdot (t - t_{lag})} - e^{-k_{i} \cdot (t - t_{lag})}) + M_{i}^{0} \cdot e^{-k_{i} \cdot t} & \text{ow.}
  \end{cases}    \label{eq4}    \tag{eq. 4} 
  \end{align}

Three additional models are derived from this general model (without degradation, without lag phase, and without degradation nor lag phase). 

Indeed, without a lag phase (i.e. :math:`t_{lag}=0`), equations 3-4 simplifies to:

.. math:: 
  \begin{align}
    X(t) = X_{0} \cdot e^{\mu \cdot t}    \label{eq5}    \tag{eq. 5} 
  \end{align}

.. math:: 
  \begin{align}
    M_{i}(t) = qM_{i} \cdot \dfrac{X_{0}}{\mu + k_{i}} \cdot (e^{\mu \cdot t} - e^{-k_{i} \cdot t}) + M_{i}^{0} \cdot e^{-k_{i} \cdot t}    \label{eq6}    \tag{eq. 6} 
  \end{align}

In the absence of degradation (i.e. :math:`k=0`), equation 4 simplifies to:

.. math:: 
  \begin{align}
    M_{i}(t) = \begin{cases}
    M_{i}^{0} & \text{if } t < t_{lag} \\
    qM_{i} \cdot \dfrac{X_{0}}{\mu} \cdot (e^{\mu \cdot (t - t_{lag})} - 1) + M_{i}^{0} & \text{ow.}
  \end{cases}    \label{eq7}    \tag{eq. 7} 
  \end{align}

In the absence of both degradation and lag (i.e. :math:`t_{lag}=0` and :math:`k=0`), equations 3-4 simplifies to:

.. math:: 
  \begin{align}
    X(t) = X_{0} \cdot e^{\mu \cdot t}    \label{eq8}    \tag{eq. 8} 
  \end{align}

.. math:: 
  \begin{align}
    M_{i}(t) = qM_{i} \cdot \dfrac{X_{0}}{\mu} \cdot (e^{\mu \cdot (t - t_{lag})} - 1) + M_{i}^{0}    \label{eq9}    \tag{eq. 9} 
  \end{align}


.. _default_dynamic_models:

Dynamic model
-----------------------

We also implemented a dynamic model where fluxes and growth are represented 
by Monod kinetics, for one substrate and one product.

Time course concentrations of 
biomass (:math:`X`), substrate (:math:`S`) and product (:math:`P`) are described by the following system of ordinary 
differential equations (ODEs):

.. math:: 
  \begin{align}
    \dfrac{dS}{dt}=-X\cdot q_{S}    \label{eq10}    \tag{eq. 10} 
  \end{align}

.. math::
  \begin{align}
    \dfrac{dX}{dt}=q_{S}\cdot yield_{biomass}    \label{eq11}    \tag{eq. 11} 
  \end{align}

.. math::
  \begin{align}
    \dfrac{dP}{dt}=q_{S}\cdot yield_{product}    \label{eq12}    \tag{eq. 12} 
  \end{align}

where :math:`q_{S}` is is the substrate uptake flux, :math:`yield_{biomass}` is the biomass yield and :math:`yield_{product}` is the product yield.

The dependence of the substrate uptake flux on the substrate concentration is expressed in this model by the Monod rate law:

.. math::
  \begin{align}
    q_{S}=q^{max}_{S}\cdot \dfrac{S}{K_{M}+S}    \label{eq13}    \tag{eq. 13} 
  \end{align}

where :math:`q^{max}_{S}` is the maximal substrate uptake rate and :math:`K_{M}` is the "half-velocity constant" (the value of :math:`S` at which :math:`\frac{\mu}{\mu_{max}}=0.5`).


User-made models
*****************


Overview
--------

Since PhysioFit 3.0.0, users can create and implement their own models to calculate fluxes and other growth parameters for any biological system. This
section explains how to write your first model, how to test the model and how to implement it
on your PhysioFit instance.

Build a template
----------------

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

This is the base template to build your model. Methods :samp:`get_params` (to initialize and return model parameters) and :samp:`simulate` (to simulate metabolite dynamics for a given set of parameters) are mandatory. Additional methods are allowed if needed (e.g. to carry out intermediary steps for the simulation).

Populate the template
---------------------

The first attribute to add in your model's :samp:`__init__` method is the model name. We strongly advise
to choose a name that helps the user understand what the model is destined to simulate. You must also add two other
attributes: the free parameters that PhysioFit will estimate & the fixed parameters provided by users. Finally, you must also call the :samp:`super().init(data)`
method to inherit the logic from the base class: ::

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

.. note:: If your model does not contain fixed parameters, you must still initialize the attribute as :samp:`None`. This is
          considered good practice.

We can now check that the model can be initialized properly. Use the block at the end of the file for
testing purposes. Here is an example of how you can test the model: ::

    if __name__ == "__main__":

        test_data = pd.DataFrame(
            {
                "time": [0, 1, 2, 3],
                "experiments": ["A", "A", "A", "A"],
                "X": [0.5, 0.8, 1.2, 1.8],
                "Glucose": [12, 11.6, 11, 10.2]
            }
        )

        model = ChildModel(data=test_data)
        print(model)

If you now run the file, you should have a standard output in your console that looks like: ::

    Selected model: Tutorial model
    Model data:
       time experiments    X  Glucose
    0     0           A  0.5     12.0
    1     1           A  0.8     11.6
    2     2           A  1.2     11.0
    3     3           A  1.8     10.2
    Experimental matrix:
    [['A' 0.5 12.0]
     ['A' 0.8 11.6]
     ['A' 1.2 11.0]
     ['A' 1.8 10.2]]
    Time vector: [0 1 2 3]
    Name vector: ['X', 'Glucose']
    Biomass & Metabolites: ['Glucose']
    Parameters to estimate: None
    Fixed parameters: None
    Bounds: None

The next step is to define the parameters (used for simulationsand optimization). PhysioFit supports two types of parameters (**parameters to estimate** and **fixed parameters**) which are detailed below.

.. _parameters_to_estimate:

Free parameters
^^^^^^^^^^^^^^^

The free parameters are the parameters that will be estimated by PhysioFit, and thus that require defining bounds and initial values
to be initialized. The list of parameters and their initial (default) values must be returned by the :samp:`get_params` method: ::

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
class of objects that handle the logic and checks. They are derived from the python :samp:`dict` base class, and as such
implement the same methods (e.g. :samp:`update`). Here is an example of how to implement the bounds: ::

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
^^^^^^^^^^^^^^^^

The fixed parameters are parameters that are known and are not estimated by PhysioFit. For example, in the case of
steady-state models that account for non enymatic degradation (see :ref:`default_steady-state_models`.), we need to provide the
degradation constant of all unstable metabolites (these constants must be measured in an independent experiment, e.g. see Peiro et al., 2019): ::

    self.fixed_parameters = {"Degradation": {
            metabolite: 2 for metabolite in self.metabolites
            }
        }

The fixed parameters must be provided as a dictionary of dictionaries, where the first level is the name of the
parameter itself (here degradation) and the second level contains the mapping of metabolite-value pairs that will be
the default values initialized (here we give a default value of 2 for every metabolite for example). Each
key of the first level is used to initialize a widget in the GUI, thus allowing users to change the corresponding 
values for the metabolites given in the second level.

Simulation function
^^^^^^^^^^^^^^^^^^^

Once the :samp:`get_params` method has been implemented, the next step is to implement the simulation function that
will be called at each iteration of the optimization process to simulate the metabolite dynamics that correspond to a 
given set of parameters (see :ref:`optimization_process` for more details).
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
    * :samp:`params_opti`: list containing the values of each parameter to estimate **in the same order as defined in the :samp:`parameters_to_estimate` dictionary** (see :ref:`parameters_to_estimate`)
    * :samp:`data_matrix`: numpy array containing the experimental data (or data with the same shape)
    * :samp:`time_vector`: numpy array containing the time points
    * :samp:`params_non_opti`: dictionary containing the fixed parameters (see :ref:`fixed_parameters`)

Now you can start writing the body of the function. For sake of clarity, we recommend unpacking parameters values from the 
list of parameters to estimate into internal variables. Th function *simulate* must return a matrix containing the simulation results, with the same shape as 
the matrix containing the experimental data. To initialize the simulated matrix, you can 
use the :samp:`empty_like` function from the numpy library: ::

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
5 and 6). See :ref:`here <testing_the_model>` for information on how to test the completed model.

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

As we can see, the function :samp:`calculate_derivative` returns the derivatives of each metabolite concentration and is used by an ODEs solver that performs the simulations. This function is thus
created within the body of the simulate function, before being called by the solver. More information on the mathematics
behind this implementation can be found :ref:`here <default_dynamic_models>`.

.. note:: The simulation function will be called a high number of times by the optimizer for parameter estimation, so optimize this function as much as possible. When possible, implement the model using analytical solution as calculations will be faster than solving numerically the corresponding ODEs.

.. _testing_the_model:

Test the model
--------------

Once you have completely populated your model file, you can now launch a round of simulations and optimizations in a
programmatic way: ::

    if __name__ == "__main__":
        from physiofit.base.io import IoHandler
        from physiofit.models.base_model import StandardDevs

        test_data = pd.DataFrame(
            {
                "time": [0, 1, 2, 3],
                "experiments": ["A", "A", "A", "A"],
                "X": [0.5, 0.8, 1.2, 1.8],
                "Glucose": [12, 11.6, 11, 10.2]
            }
        )

        io = IoHandler()
        model = ChildModel(data=test_data)
        model.get_params()
        fitter = io.initialize_fitter(
            model.data,
            model=model,
            mc=True,
            iterations=100,
            sd=StandardDevs({"X": 1, "Glucose": 1}),
            debug_mode=True
        )
        fitter.optimize()
        fitter.monte_carlo_analysis()
        fitter.khi2_test()
        print(fitter.parameter_stats)

This will return the calculated flux values and associated statistics.

To test the integration of the model into the GUI, copy the :file:`.py` file
in the folder :file:`models` of PhysioFit directory. You can get the path towards this folder by opening a python
kernel in your dedicated environment and initializing an IoHandler ::

    from physiofit.base.io import IoHandler
    io_handler = IoHandler()
    print(io_handler.get_local_model_folder())

.. note:: The model file name must follow the naming convention :file:`model_[model number].py`. If the last model in the list
          is the :file:`model_5.py`, the next one should be named :file:`model_6.py`.

You can now launch PhysioFit's GUI, load a data file corresponding to the new model, select the model, and run flux calculation. In case of errors, 
have a look to the error message and correct the code.

.. note:: We would be happy to broaden the types of models shipped with PhysioFit. If you have developed a new model, it might be 
          usefull and valuable to the fluxomics community! Please, keep in touch with us to discuss on the model and see if we can include your 
          model in the built-in models shipped with PhysioFit! :)