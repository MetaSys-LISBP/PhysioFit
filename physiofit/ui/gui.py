import json
import tkinter as tk
from ast import literal_eval
from copy import copy
from pathlib import Path
from threading import Thread
from tkinter import filedialog

import pandas as pd
import requests
import streamlit as st

import physiofit
from physiofit.base.io import IoHandler, DEFAULTS


class App:
    """
    Physiofit Graphical User Interface
    """

    def __init__(self):

        self.defaults = copy(DEFAULTS)
        self.select_menu = None
        self.io_handler = None
        self.update_info = None

    def start_app(self):
        """Launch the application"""
        st.title("Welcome to PhysioFit")
        self.update_info = st.empty()
        # thread = Thread(target=self.get_last_version)
        # thread.start()
        self.check_uptodate()
        self.select_menu = st.selectbox(
            "Select a task to execute",
            ["Calculate extracellular fluxes", "Calculate degradation constant"]
        )
        if self.select_menu == "Calculate extracellular fluxes":
            self.io_handler = IoHandler()
            self._build_flux_menu()
        else:
            st.header("Implementation in progress...")

    def check_uptodate(self, package_name="physiofit"):
        """Compare installed and most recent Physiofit versions."""

        try:
            response = requests.get(f"https://pypi.org/pypi/{package_name}/json")
            latest_version = response.json()['info']['version']
        except Exception:
            latest_version = False
        if not latest_version:
            st.info(
                f"Error while retrieving latest version. Please check manually. "
                f"Current version: {physiofit.__version__}"
            )
        elif latest_version != physiofit.__version__:
            st.info(
                f'New version available ({latest_version}). \n You can update PhysioFit with: "pip install --upgrade '
                f'physiofit" \n Check the documentation for more information.'
            )

    def _build_flux_menu(self):
        """Build the starting menu with the data upload button"""

        self.data_file = st.file_uploader(
            "Load a data file or a json configuration file",
            key="data_uploader",
            accept_multiple_files=False
        )

        if self.data_file:
            self._initialize_opt_menu()

    def _initialize_opt_menu(self):
        """
        Initialize all the options. If a json file is given as input, options are parsed from it. If the input is
        tsv, the defaults are used instead.
        """

        file_extension = self.data_file.name.split(".")[1]
        input_values = self.defaults
        if file_extension == "json":
            config = IoHandler.read_json_config(self.data_file)
            input_values.update(config)
        elif file_extension != "tsv":
            raise KeyError(
                f"Wrong input file format. Accepted formats are tsv for data files or json for configuration "
                f"files. Detected file: {self.data_file.name}")
        else:
            data = pd.read_csv(self.data_file, sep="\t")
            try:
                self.defaults["weight"].update({"X": 0.2})
                for col in data.columns[2:]:
                    self.defaults["weight"].update({col: 0.5})
            except Exception:
                raise

        submitted = self._initialize_opt_menu_widgets(input_values, file_extension)

        if submitted:
            if file_extension == "tsv":
                self.io_handler.data = data
                self.io_handler.data = self.io_handler.data.sort_values("time", ignore_index=True)
                self.io_handler.names = self.io_handler.data.columns[1:].to_list()
                kwargs = self._build_fitter_kwargs()
                self.io_handler.initialize_fitter(kwargs)
            if file_extension == "json":
                st.session_state.submitted = True
                final_json = self._build_internal_json(input_values["path_to_data"])
                self.io_handler.launch_from_json(final_json)
            self.io_handler.fitter.optimize()
            if self.mc:
                self.io_handler.fitter.monte_carlo_analysis()
            self.io_handler.fitter.khi2_test()
            outputs = ["data", "plot", "pdf"]
            self.io_handler.local_out(*outputs)
            st.write(f"Run is finished. Check {self.io_handler.res_path} for the results.")

    def _initialize_opt_menu_widgets(self, input_values, file_extension):

        expand_basic_settings = st.expander("Basic settings", expanded=True)
        with expand_basic_settings:
            self.t_lag = st.checkbox(
                "Lag",
                key="lag_check"
            )
            self.deg_check = st.checkbox(
                "Degradation",
                key="deg_check"
            )
            enable_deg = False if self.deg_check else True
            self.deg = st.text_input(
                "Degradation constants",
                value={},
                help="Dictionary containing the (first-order) degradation constant of each metabolite. "
                     "Format: {'met1': value, 'met2': value, ... 'metn': value}",
                disabled=enable_deg
            )
            self.mc = st.checkbox(
                "Sensitivity analysis (Monte Carlo)",
                value=True,
                help="Determine the precision on estimated fluxes by Monte Carlo sensitivity analysis."
            )
            enable_mc = False if self.mc else True
            self.iterations = st.number_input(
                "Number of iterations",
                value=input_values["iterations"],
                help="Number of iterations for the Monte Carlo analysis.",
                disabled=enable_mc
            )

            if file_extension == "tsv":

                # Set up tkinter for directory chooser
                root = tk.Tk()
                root.withdraw()

                # Make folder picker dialog appear on top of other windows
                root.wm_attributes('-topmost', 1)

                # Initialize folder picker button and add logic
                clicked = st.button("Select output data directory", key="clicker")
                if clicked:

                    # Initialize home path from directory selector and add to session state
                    st.session_state.home_path = Path(st.text_input(
                        "Selected output data directory:", filedialog.askdirectory(master=root)
                    ))
                    if st.session_state.home_path == Path(".") or not st.session_state.home_path.exists():
                        raise RuntimeError("Please provide a valid path")
                    self.io_handler.home_path = copy(st.session_state.home_path)

                elif hasattr(st.session_state, "home_path"):

                    self.io_handler.home_path = Path(st.text_input(
                        "Selected output data directory:", st.session_state.home_path
                    ))
                    if self.io_handler.home_path == Path(".") or not self.io_handler.home_path.exists():
                        raise RuntimeError("Please provide a valid path")

                    # Initialize the result export directory
                    self.io_handler.res_path = self.io_handler.home_path / (self.io_handler.home_path.name + "_res")
                    if not clicked:
                        if not self.io_handler.res_path.is_dir():
                            self.io_handler.res_path.mkdir()

        # Build the form for advanced parameters
        form = st.form("Parameter_form")
        with form:
            expand_run_params = st.expander("Advanced settings")
            with expand_run_params:
                self.vini = st.text_input(
                    "Initial flux values",
                    value=input_values["vini"],
                    key="vini",
                    help="Select an initial value of fluxes to estimate. Default: 0.2"
                )
                self.weight = st.text_input(
                    "Standard deviation on measurements",
                    value=input_values["weight"],
                    help="Standard deviation on the measurements. If empty, default is 0.02 for biomass and"
                         " 0.5 for metabolites"
                )
                self.conc_met_bounds = st.text_input(
                    "Bounds on initial metabolite concentrations",
                    value=input_values["conc_met_bounds"],
                    help="Bounds for the initial concentrations of the metabolites (Mi0). "
                         "These values correspond to the lowest and highest initial concentration of metabolites, "
                         "this range should include the actual values. Defaults: [1e-06, 50]"
                )
                self.flux_met_bounds = st.text_input(
                    "Bounds on fluxes",
                    value=input_values["flux_met_bounds"],
                    help="Bounds for metabolic fluxes (qM). "
                         "These values correspond to the lowest and highest fluxes, this range should include the "
                         "actual value. Defaults: [0.01, 50]"
                )
                self.conc_biom_bounds = st.text_input(
                    "Bounds on initial biomass concentration",
                    value=input_values["conc_biom_bounds"],
                    help="Bounds for initial concentrations of the biomass (X0). "
                         "These values correspond to the lowest and highest (initial) biomass concentration, this "
                         "range should include the actual value. Defaults: [1e-06, 50]"
                )
                self.flux_biom_bounds = st.text_input(
                    "Bounds on growth rate",
                    value=input_values["flux_biom_bounds"],
                    help="Bounds for growth rate (µ). "
                         "These values correspond to the lowest and highest growth rates, this range should include "
                         "the actual value. Defaults: [0.01, 2]"
                )
                self.debug_mode = st.checkbox(
                    "Verbose logs",
                    help="Useful in case of trouble. Join it to the issue on github."
                )
            submitted = st.form_submit_button("Run flux calculation")
        return submitted

    def _build_fitter_kwargs(self):

        kwargs = {
            "vini": float(self.vini),
            "weight": literal_eval(self.weight),
            "conc_met_bounds": tuple(literal_eval(self.conc_met_bounds)),
            "flux_met_bounds": tuple(literal_eval(self.flux_met_bounds)),
            "conc_biom_bounds": tuple(literal_eval(self.conc_biom_bounds)),
            "flux_biom_bounds": tuple(literal_eval(self.flux_biom_bounds)),
            "t_lag": self.t_lag,
            "deg": literal_eval(self.deg),
            "mc": self.mc,
            "iterations": self.iterations,
            "debug_mode": self.debug_mode,
        }
        return kwargs

    def _build_internal_json(self, path_to_data):

        final_json = json.dumps({
            "vini": float(self.vini),
            "weight": literal_eval(self.weight),
            "conc_met_bounds": literal_eval(self.conc_met_bounds),
            "flux_met_bounds": literal_eval(self.flux_met_bounds),
            "conc_biom_bounds": literal_eval(self.conc_biom_bounds),
            "flux_biom_bounds": literal_eval(self.flux_biom_bounds),
            "t_lag": self.t_lag,
            "deg": literal_eval(self.deg),
            "mc": self.mc,
            "iterations": self.iterations,
            "debug_mode": self.debug_mode,
            "path_to_data": path_to_data
        })
        return final_json


if __name__ == "__main__":
    physiofit_app = App()
    physiofit_app.start_app()
