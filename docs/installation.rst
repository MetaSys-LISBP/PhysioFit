Installation
============

Installation
-----------------

PhysioFit requires Python 3.6 or higher and can run on most systems supporting Python3 (Windows, MacOS and Linux). If you do not have a Python environment
configured on your computer, we recommend that you follow the instructions
from `Anaconda <https://www.anaconda.com/download/>`_.

To install PhysioFit using Python's built-in installer, you can just run the following command in a terminal:

.. code-block:: bash

    pip install physiofit

.. tip::  We recommend the creation of isolated environments for each python tool you install in your system using the python built-in `venv <https://docs.python.org/3/library/venv.html>`_ package or `Anaconda <https://www.anaconda.com/products/individual>`_.

If this method does not work, you should ask your local system administrator or
the IT department "how to install a Python 3 package from PyPi" on your computer.

PhysioFit is freely available and is distributed under open-source license at http://github.com/MetaSys-LISBP/ .


Alternatives & Updates
----------------------

If you know that you do not have permission to install software system-wide, you can install PhysioFit into your user directory using the `--user` flag:

.. code-block::

    pip install --user physiofit

This does not require any special privileges.

Once the package is installed, you can update it using the following command:

.. code-block::

    pip install -U physiofit

Alternatively, you can also download all sources in a tarball from `GitHub <https://github.com/MetaSys-LISBP/PhysioFit>`_, but it will be more difficult to update PhysioFit later on.
