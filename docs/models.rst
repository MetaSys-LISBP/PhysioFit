Models
=========

Default models
***************

Overview
-----------------

Initial concentrations of species (i.e. biomass - X - and metabolites - M_i - at t=0) and fluxes
(exchange fluxes - q_(M_i ) - and growth rate - µ -) are estimated by fitting time-course measurements of
metabolite and biomass concentrations. Simulations, i.e. calculation of X(t) and M_i(t), are performed using analytical
functions, as detailed below in the Model section.

Concentrations of species (biomass and metabolites) at t=0 and exchange fluxes are estimated by fitting time-course
measurements of metabolite and biomass concentrations. Different models are shipped by default with PhysioFit. This
section of the documentation details each of these models.

Flux values given by PhysioFit correspond the best fit. A global sensitivity analysis (Monte-Carlo approach) is
available to evaluate the precision of the estimated fluxes (mean, median, standard deviation, 95% confidence
intervals), plots are generated for visual inspection of the fitting quality, and a χ² test is performed to assess the
statistical goodness of fit.

.. _steady-state-models:

Steady-state models
--------------------

The steady state models implemented in PhysioFit may account for i) non enzymatic degradation of some metabolites and
ii) growth lag. These models can be described using the following system of ordinary differential equations:

.. image:: _static/equations/eq1.png

.. image:: _static/equations/eq2.png

with qMi being positive (negative) when M_i is produced (consumed). The sign of q_(M_i ) can thus be used to
automatically identify products and substrates in high throughput workflows for automated functional analysis of
metabolic systems.

Integrating equations 1-2 provides the following analytical functions:

.. image:: _static/equations/eq3.png

.. image:: _static/equations/eq4.png

In the absence of a lag phase (i.e. t_lag=0), equations 3-4 simplifies to:

.. image:: _static/equations/eq5.png

.. image:: _static/equations/eq6.png

In the absence of degradation (i.e. k = 0), eq. 4 simplifies to:

.. image:: _static/equations/eq7.png

In the absence of both degradation and lag (i.e. t_lag=0 and k=0), equations 3-4 simplifies to:

.. image:: _static/equations/eq8.png

.. image:: _static/equations/eq9.png


Dynamic models
----------------

To be implemented...

Flux calculation
-----------------

A model containing all metabolites is constructed automatically by PhysioFit from the input data, as detailed above.
Model parameters (i.e. fluxes, growth rate, and initial concentrations of biomass and metabolites) are estimated using
Scipy's Differential evolution method, and the best solution is polished using the L-BFGS-B method (see
`scipy.optimize <https://docs.scipy.org/doc/scipy/reference/optimize.html>`_ for more
information on the optimization method), by minimizing the following cost function:

.. image:: _static/equations/eq10.png

where *sim* is the simulated data, *meas* denotes measurements, and *weight* is a weighting factor (typically, one standard
deviation on measurements).


User-made models
*****************

Overview
---------

Since PhysioFit 3.0.0, it is possible to create and implement models into the tool using your own python3 code. This
section of documentation will explain how to write your first model, how to test your model and implement it
into your local installation. A section will also describe how to submit your code for integration onto the
`Workflow4Metabolomics <https://workflow4metabolomics.usegalaxy.fr/>`_ platform for use in fluxomics workflows for example.

# TODO: Add method section for flux calculation & sensitivity analysis and have it be first
# TODO: Change init growth rate to BM_0
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
                "init_growth_rate": 1,
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
                "init_growth_rate": 1,
                "growth_rate": 1
            }

            # Instantiate the bounds object

            self.bounds = Bounds(
                {
                    "init_growth_rate": (1e-3, 10),
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


Fixed parameters
""""""""""""""""

The fixed parameters are parameters that are given as constants in the model equations. For example, in the case of
steady-state models that account for non enymatic degradation (see :ref:`steady-state-models`.), we need to give
the unstable metabolite a constant that will define it's rate of degradation ::

    self.fixed_parameters = {"Degradation": {
            metabolite: 2 for metabolite in self.metabolites
            }
        }

The different fixed parameters are given in a dictionnary of dictionnaries, where the first level is the name of the
parameter itself (here degradation) and the second level contains the mapping of metabolite-value pairs that will be
the default values initialized on launch(here we give a default value of 0 for every metabolite for example).


