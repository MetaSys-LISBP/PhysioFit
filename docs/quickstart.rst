Quick start
==========

In this section we will explain how to launch your first job once PhysioFit has been installed onto your system.

.. seealso::
   * If you have already used PhysioFit and are looking for a more in-depth tutorial, check out the :doc:`usage` section.

   * The mathematical formalism implemented in PhysioFit to calculate fluxes can be found in the :doc:`method` section.

Graphical user interface
--------------------------------------

To launch the Graphical User Interface, type in a terminal (e.g. Anaconda Prompt if installed on Windows):

.. code-block:: bash

  physiofit
 
If you installed the package in a specific environment, make sure to activate this environment before starting PhysioFit.

PhysioFit interface will open in a new browser window.

.. image:: _static/interface.jpg
   :scale: 85%

Select an input 
file (which can be a .tsv file containing the data or a json configuration file containing the run parameters and
a path towards the data, see :doc:`usage` for more details), modify the calculation parameters according to your experiment, 
and click on :samp:`Run flux calculation`. PhysioFit proceeds automatically to the flux calculation and display its progress and important messages. The output of the calculations (i.e. fluxes and associated confidence intervals) will be written in a text file along a log file.

Command line interface
----------------------

Implementation in progress...

Library
-------

PhysioFit is also available as a library (a Python module) that you can import directly in your Python
scripts:

.. code-block:: python

  import physiofit

.. seealso::  Have a look at our :ref:`library showcase <Library documentation>` if you are interested into this experimental feature.