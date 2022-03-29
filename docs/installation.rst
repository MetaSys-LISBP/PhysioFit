Installation
============

Installation
-----------------

PhysioFit requires Python 3.6 or higher and can run on most systems supporting Python3 (Windows, MacOS and Linux). If you do not have a Python environment
configured on your computer, we recommend that you follow the instructions
from `Anaconda <https://www.anaconda.com/download/>`_.

To install PhysioFit 2.0 using Python's built-in installer, you can just run the following command in a terminal:

.. code-block:: bash

    pip install physiofit

We recommend you create isolated environments for each python tool you install in your system using the python built-in
`venv <https://docs.python.org/3/library/venv.html>`_ package or `Anaconda <https://www.anaconda.com/products/individual>`_.

If this method does not work, you should ask your local system administrator or
the IT department "how to install a Python 3 package from PyPi" on your computer.

PhysioFit is freely available and is distributed under open-source license at http://github.com/MetaSys-LISBP/ .


Alternatives & Updates
----------------------

To install PhysioFit 2.0 for the current user account only, run the following command:

.. code-block::

    pip install --user physiofit

This will ensure that the package and dependencies are installed to your user's home directory instead of a system
directory, which does not require any special privileges.

Once the package is installed, you can update it using the following command:

.. code-block::

    pip install -U physiofit

