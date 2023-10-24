Quick start
============

In this section we will explain how to launch your first job once PhysioFit has been installed onto your system.

.. seealso::
   * If you have already used PhysioFit and are looking for a more in-depth tutorial, check out the :doc:`usage` section.

Graphical user interface
--------------------------------------

To open the Graphical User Interface, type in a terminal (e.g. Anaconda Prompt if installed on Windows):

.. code-block:: bash

  physiofit
 
If you have installed the package in a specific environment, make sure to activate this environment before starting PhysioFit.

PhysioFit interface will open in a new browser window:

.. image:: _static/interface.jpg
   :scale: 60%

Select an input file (which can be a :file:`tsv` file containing the data or a :file:`yaml` configuration file containing the run
parameters and a path towards the data, see :doc:`usage` for more details), select a model, modify the calculation parameters according
to your data, and click on :samp:`Run flux calculation`. PhysioFit proceeds automatically to the flux calculation
and display its progress and possibly important messages such as errors. The output of the calculations (i.e. fluxes and associated statistics) 
will be written in a text file as will the statistical test results, while plots will be generated individually for each metabolite (:file:`svg` files) and combined in a 
multi-page :file:`pdf` file. If multiple experiments were included in the input data, a summary (:file:`csv` file) 
will also be generated. See :ref:`outputs_ref` for more details.

Command line interface
----------------------

To process your data, type in a terminal:

.. code-block:: bash

  physiofit [command line options]

Here after the available options with their full names are enumerated and detailed.

.. argparse::
   :module: physiofit.ui.cli
   :func: args_parse
   :prog: physiofit
   :nodescription:

Library
-------

PhysioFit is also available as a library (a Python module) that you can import directly in your Python
scripts:

.. code-block:: python

  import physiofit

.. seealso::  Have a look at our :ref:`API <Library documentation>` if you are interested in this feature.